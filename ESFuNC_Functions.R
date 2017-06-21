###
# Ensemble of Segmented Functional Nonparametric Classifiers
# Robert Buscaglia, Nichola C. Garbett, Yiannis Kamarianakis
# June 20, 2017
#
# Functions and Packages
# File contains all required functions and packages for running the
# ESFuNC algorithm.  Each function has been given a brief description.
# Please see ESFuNC_Analysis_Template.R for details on performing the
# algorithm.

### Required Packages
# Checks for and installs missing packages then loads all packages.

list.of.packages <- c("fda.usc", "foreach", "doParallel", "abind")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) 
  {
    cat("Missing Packages are Installing... \n")
    install.packages(new.packages)
  }

library("fda.usc")
library("foreach")
library("doParallel")
library("abind")

### Kernels
# Provide weights for nonparametric classifiers. Required for probability calculations.  
# Must be written so that they can work with both KNN and Parzen Window.  Any additional 
# kernels can be prepared but should retain this form to ensure functionality.           

kern.tri<-function(x, h=1, method="wknn")
{
  if(method=="wknn")
  {
    temp.1<-x-min(x)
    temp.2<-temp.1/max(temp.1)
    result<-2*(1-temp.2)
  }
  
  if(method=="kernel")
  {
    result<-ifelse(x<=h, (2/h)*(1-(x/h)), 0)
  }
  
  return(result)
}


kern.norm<-function(x, h=1, method="wknn")
{
  if(method=="wknn")
  {
    temp.1<-x-min(x)
    result<-2*dnorm(temp.1)
  }
  
  if(method=="kernel")
  {
    result<-2*dnorm(x, 0, h)
  }
  
  return(result)
}

kern.unif<-function(x, h=1, method="wknn")
{
  if(method=="wknn")
  {
    result<-dunif(x, min = min(x), max=max(x))
  }
  
  if(method=="kernel")
  {
    result<-ifelse(x<=h, 1/h, 0)
  }
  
  return(result)
}

### createFolds
# This is the createFolds function from the package 'caret'.
# This is  used to create stratified folds for cross validation.

createFolds<-function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}


### fold.creator
# A default setting for creating multiple sets of folds if not given manually.

fold.creator<-function(classes, folds, trials)
{
  folds.list<-list()
  for(j in  1:trials)
    folds.list[[j]]<-createFolds(factor(classes), k=folds, list=FALSE)
  return(folds.list)
}

### create.cv.sum
# A function for aiding in the visualization and summarization of final cross validation results. Takes
# the results of running any of the cross validation functions and produces summary statistics.
# cv.to.plot should be at least one set of cross validation results.  This is a common output from
# cross validation functions.

create.cv.sum<-function(cv.to.plot)
{
  n.size<-nrow(as.matrix(cv.to.plot))
  med.temp<-apply(cv.to.plot, 2, median)
  mad.temp<-apply(cv.to.plot, 2, mad)
  mean.temp<-apply(cv.to.plot, 2, mean)
  sd.temp<-apply(cv.to.plot, 2, sd)
  se.temp<-sd.temp/sqrt(n.size)
  ul.temp<-mean.temp+2*se.temp
  ll.temp<-mean.temp-2*se.temp
  return(rbind(med.temp, mad.temp, ul.temp, mean.temp, sd.temp, ll.temp))
}

### create.segments
# Creates unique fdata objects for each segment determined by (total.segments)
# total.segments = 1 refers to the undisturbed fdata object.
# fdata.object should be fdata object being analyzed (any derivative)
# total.segments : how many partitions are to be made
# indent : discrete point indentations based on argvals.  used to remove numerical noise.
# density : grid mesh point density to control integration accuracy.

create.segments<-function(fdata.object, total.segments, indent, density=NULL)
{
  argvals=fdata.object$argvals
  size<-length(argvals)
  start<-argvals[1+indent]
  end<-argvals[size-indent]
  
  ifelse(is.null(density), int.density<-size-2*indent, int.density<-density)
  temp<-0:total.segments
  seg.dist<-(end-start)/total.segments
  lims.temp<-start+temp*seg.dist
  
  segs.list<-list()
  for(j in 1:total.segments)
  {
    seq.temp<-seq(lims.temp[j], lims.temp[j+1], (lims.temp[j+1]-lims.temp[j])/(int.density-1))
    fd.temp<-eval.fd(seq.temp, fdata2fd(fdata.object))
    segs.list[[j]]<-fdata(t(fd.temp), argvals=seq.temp)
  }
  
  return(segs.list)
}

### sort.models
# Sorts grid search results and orders results such that segments are minimized and neighbors are maximized
# grid.output : resulting object from running sequential grid search
# models.returned : how many models to be viewed as output
# do.print : control whether resulting ordered table of values is printed to screen.  Used in other functions.

sort.models<-function(grid.output, models.returned, do.print=TRUE)
{
  size<-dim(grid.output)[1]*dim(grid.output)[2]
  best.models<-arrayInd(order(grid.output, decreasing = TRUE), dim(grid.output))
  best.models.mat<-matrix(nrow=models.returned, ncol=3)
  best.models.mat.temp<-matrix(nrow=size, ncol=3)
  colnames(best.models.mat)<-c("Neighbors", "Segments", "Accuracy")
  for(j in 1:size)
  {
    best.models.mat.temp[j,]<-c(rownames(grid.output)[best.models[j,1]],colnames(grid.output)[best.models[j,2]],as.numeric(round(grid.output[best.models[j,1],best.models[j,2]],6)))
  }
  
  acc.levels<-as.numeric(levels(factor(best.models.mat.temp[,3])))
  for(j in 1:length(acc.levels))
  {
    set.1<-which(best.models.mat.temp[,3]==acc.levels[j])
    fac.levels<-as.numeric(levels(factor(best.models[set.1,2])))
    fac.len<-length(fac.levels)
    for(k in 1:fac.len)
    {
      temp.1<-which(best.models[set.1,2]==fac.levels[k])
      best.models[set.1[temp.1]]=sort(best.models[set.1[temp.1]], decreasing = TRUE)
    }
  }
  
  for(j in 1:models.returned)
  {
    best.models.mat[j,]<-c(rownames(grid.output)[best.models[j,1]],colnames(grid.output)[best.models[j,2]],as.numeric(round(grid.output[best.models[j,1],best.models[j,2]],6)))
  }
  if(do.print) print(best.models.mat, quote=FALSE)
  if(!do.print) return(best.models.mat)
}

### Stepwise Procedures ###

### all.segs.ensemble
# Does not do any stepwise selection.  Uses all segments in final model.
# step.array : array of probabilities determined from distances, and dependent on which method and kernel chosen.
#              Calculated frequently in other functions.
# segment.asccuracies : accuracy of each segment included in step.array.  Used for weighted combinations.
# classes : classification identifiers for set being analyzed.                                             
# seg.weight : Set to TRUE if combining segment probabilities should be weighted by individual segment accuracy
# thresh : stepwise accuracy improvement threshold.  Does not affect all.segs.ensemble.                     
# do.par : Set to TRUE to run calculations in parallel.                                                    
# cores : number of cores to use during parallel calculations.                                            

all.segs.ensemble<-function(step.array, segment.accuracies, classes, seg.weight=FALSE, thresh=0.0001, do.par=FALSE, cores=2)
{
  n.objs<-dim(step.array)[1]
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  knn.n<-dim(step.array)[3]
  segment.n<-dim(step.array)[4]
  seg.seq<-seq(1,segment.n,1)
  # Add Check that Dimensions are correct
  
  ens.acc<-numeric()
  ens.segs.used<-list()
  probs.ens<-array(dim=c(n.objs, length(class.levels.index), knn.n))
  classes.ens<-matrix(nrow=n.objs, ncol=knn.n)
  
  for (k in 1:knn.n)  # Will be changed to bandwith
  {
    segs.used<-seg.seq
    
    probs.update<-matrix(nrow=n.objs, ncol=length(class.levels.index))
    for (q in class.levels.index)
    {
      ifelse(seg.weight, 
             probs.update[,q]<-rowSums(t(segment.accuracies[k,segs.used]*t(as.matrix(step.array[,q,k,segs.used])))/sum(segment.accuracies[k,segs.used])), 
             probs.update[,q]<-rowSums(as.matrix(step.array[,q,k,segs.used]))/length(segs.used)
      )
    }
    
    #print(probs.update)
    
    est.classes.update<-numeric()
    for (q in 1:n.objs)
    {
      max.classes<-class.levels[as.numeric(which(probs.update[q,]==max(probs.update[q,])))]
      est.classes.update[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))      
    }
    
    acc.new<-mean(est.classes.update==classes)
    
    ens.acc[k]<-acc.new
    ens.segs.used[[k]]<-segs.used
    probs.ens[,,k]<-probs.update
    classes.ens[,k]<-est.classes.update
  }
  return(list(ens.accuracies=ens.acc, ens.segments=ens.segs.used, ens.probs=probs.ens, ens.classes=classes.ens))
}

### forward.ensemble
# Evaluates the sequential addition of segments to the ensemble, evaluating if inclusion of a new segment
# improves leave-one-out accuracy.
# step.array : array of probabilities determined from distances, and dependent on which method and kernel chosen.
#              Calculated frequently in other functions.
# segment.asccuracies : accuracy of each segment included in step.array.  Used for weighted combinations.
# classes : classification identifies for set being analyzed.
# seg.weight : Set to TRUE if combining segment probabilities should be weighted by individual segment accuracy.
# thresh : stepwise accuracy improvement threshold.
# do.par : Set to TRUE to run calculations in parallel.
# cores : number of cores to use during parallel calculations.

forward.ensemble<-function(step.array, segment.accuracies, classes, seg.weight=FALSE, thresh=0.0001, do.par=FALSE, cores=2)
{
  #print("forward")
  #browser()
  n.objs<-dim(step.array)[1]
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  knn.n<-dim(step.array)[3]
  segment.n<-dim(step.array)[4]  #changed this!!!
  seg.seq<-seq(1,segment.n,1)
  # Add Check that Dimensions are correct
  
  ens.acc<-numeric()
  ens.segs.used<-list()
  probs.ens<-array(dim=c(n.objs, length(class.levels.index), knn.n))
  classes.ens<-matrix(nrow=n.objs, ncol=knn.n)
  
  if(do.par)
  {
    use.cores<-min(knn.n, max.cores)
    cl.temp<-makeCluster(use.cores)
    registerDoParallel(cl.temp)
    k.result<-foreach(k=1:knn.n, .combine=append) %dopar% #this can be parallelized
    {
      acc.old<-0 # starting value for while
      
      max.vals<-which(segment.accuracies[k,]==max(segment.accuracies[k,]))
      max.seg<-ifelse(length(max.vals)==1, max.vals, sample(max.vals,1))
      acc.new<-segment.accuracies[k,max.seg]
      segs.used<-max.seg
      segs.left<-seg.seq[-segs.used]
      
      probs.update<-step.array[,,k,segs.used]
      est.classes.update<-numeric()
      for(q in 1:n.objs)
      {
        max.classes<-class.levels[as.numeric(which(probs.update[q,]==max(probs.update[q,])))]
        est.classes.update[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
      }
      
      
      # cylce through remaining segs
      
      while((acc.new > acc.old+thresh) & (length(segs.left)!=0))
      {
        acc.old<-acc.new
        
        acc.update<-0
        seg.update<-NULL
        est.classes.update<-NULL
        
        
        for(j in segs.left)
        {
          
          probs.temp<-matrix(nrow=n.objs, ncol=length(class.levels.index))
          for(q in class.levels.index) 
          {
            ifelse(seg.weight,
                   probs.temp[,q]<-rowSums(t(segment.accuracies[k,c(segs.used,j)]*t(as.matrix(step.array[,q,k,c(segs.used,j)])))/sum(segment.accuracies[k,c(segs.used,j)])),
                   probs.temp[,q]<-rowSums(as.matrix(step.array[,q,k,c(segs.used,j)]))/length(c(segs.used,j))
            )
            if(anyNA(probs.temp[,q])) probs.temp[,q]<-rep(0, n.objs)
          }
          
          
          est.classes.temp<-numeric()
          for(q in 1:n.objs)
          {
            max.classes<-class.levels[as.numeric(which(probs.temp[q,]==max(probs.temp[q,])))]
            est.classes.temp[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          
          acc.temp<-mean(est.classes.temp==classes)
          if(acc.temp>acc.update)
          {
            acc.update<-acc.temp
            seg.update<-j
            est.classes.update<-est.classes.temp
            probs.update<-probs.temp
          }
        }
        
        
        if(acc.update>acc.old+thresh)
        {
          acc.new<-acc.update
          segs.used<-c(segs.used,seg.update)
          segs.left<-seg.seq[-segs.used]
        }
      }
      
      return(list(list(ens.acc.out=acc.new, ens.segs.used.out=segs.used, probs.ens.out=probs.update, classes.ens.out=est.classes.update)))
    }
    stopCluster(cl.temp)
  }
  
  if(!do.par)
  {
    k.result<-foreach(k=1:knn.n, .combine=append) %do% #this can be parallelized
    {
      acc.old<-0 # starting value for while
      
      max.vals<-which(segment.accuracies[k,]==max(segment.accuracies[k,]))
      max.seg<-ifelse(length(max.vals)==1, max.vals, sample(max.vals,1))
      acc.new<-segment.accuracies[k,max.seg]
      segs.used<-max.seg
      segs.left<-seg.seq[-segs.used]
      
      probs.update<-step.array[,,k,segs.used]
      est.classes.update<-numeric()
      for(q in 1:n.objs)
      {
        max.classes<-class.levels[as.numeric(which(probs.update[q,]==max(probs.update[q,])))]
        est.classes.update[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
      }
      
      
      # cylce through remaining segs
      
      while((acc.new > acc.old+thresh) & (length(segs.left)!=0))
      {
        acc.old<-acc.new
        
        acc.update<-0
        seg.update<-NULL
        est.classes.update<-NULL
        
        
        for(j in segs.left)
        {
          
          probs.temp<-matrix(nrow=n.objs, ncol=length(class.levels.index))
          for(q in class.levels.index) 
          {
            ifelse(seg.weight,
                   probs.temp[,q]<-rowSums(t(segment.accuracies[k,c(segs.used,j)]*t(as.matrix(step.array[,q,k,c(segs.used,j)])))/sum(segment.accuracies[k,c(segs.used,j)])),
                   probs.temp[,q]<-rowSums(as.matrix(step.array[,q,k,c(segs.used,j)]))/length(c(segs.used,j))
            )
            if(anyNA(probs.temp[,q])) probs.temp[,q]<-rep(0, n.objs)
          }
          
          
          est.classes.temp<-numeric()
          for(q in 1:n.objs)
          {
            max.classes<-class.levels[as.numeric(which(probs.temp[q,]==max(probs.temp[q,])))]
            est.classes.temp[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          
          acc.temp<-mean(est.classes.temp==classes)
          if(acc.temp>acc.update)
          {
            acc.update<-acc.temp
            seg.update<-j
            est.classes.update<-est.classes.temp
            probs.update<-probs.temp
          }
        }
        
        
        if(acc.update>acc.old+thresh)
        {
          acc.new<-acc.update
          segs.used<-c(segs.used,seg.update)
          segs.left<-seg.seq[-segs.used]
        }
      }
      
      return(list(list(ens.acc.out=acc.new, ens.segs.used.out=segs.used, probs.ens.out=probs.update, classes.ens.out=est.classes.update)))
    }
  }
  
  for(k in 1:knn.n)
  {
    ens.acc[k]<-k.result[[k]]$ens.acc.out
    ens.segs.used[[k]]<-k.result[[k]]$ens.segs.used.out
    probs.ens[,,k]<-k.result[[k]]$probs.ens.out
    classes.ens[,k]<-k.result[[k]]$classes.ens.out
  }
  
  return(list(ens.accuracies=ens.acc, ens.segments=ens.segs.used, ens.probs=probs.ens, ens.classes=classes.ens))
}

### calc.distance.array
# Creates an array of distances for an fdata object for a provided number of segments.
# Segments are first constructed using indent and density inputs.
# Distances are calculated between each observation in the fdata object for each segment. 
# Can be done (and should be done) in parallel.
# fdata.object : data set to be analyzed.
# indent : discrete point indentations to remove numerical derivative noise.
# density : integration grid mesh to control accuracy.
# total.segments : number of partitions, where 1 defaults to unparitioned data.
# do.par : Set to TRUE to run calculations in parallel.
# cores : number of cores to use during parallel calculations.

calc.distance.array<-function(fdata.object, indent=0, density=NULL, total.segments=1, do.par=FALSE, max.cores=2)
{
  argvals<-fdata.object$argvals
  segs<-create.segments(fdata.object, total.segments, indent, density)
  
  acomb <- function(...) abind(..., along=3)
  
  if(do.par)
  {
    use.cores<-min(total.segments, max.cores)
    cl.temp<-makeCluster(use.cores)
    registerDoParallel(cl.temp)
    seg.distance.out<-foreach(q=1:total.segments, .packages = "fda.usc", .combine='acomb', .multicombine = TRUE) %dopar%
    {
      x<-metric.lp(segs[[q]])
      return(x)
    }
    stopCluster(cl.temp)
  }
  
  if(!do.par)
  {
    seg.distance.out<-foreach(q=1:total.segments, .packages = "fda.usc", .combine='acomb', .multicombine = TRUE) %dopar%
    {
      x<-metric.lp(segs[[q]])
      return(x)
    }
  }
  
  dim.1<-dim(seg.distance.out)[1]
  dim.2<-dim(seg.distance.out)[2]
  dim.3<-total.segments
  export.array<-array(dim=c(dim.1, dim.2, dim.3))
  export.array[,,(1:total.segments)]<-seg.distance.out
  
  return(export.array)
}

### segment.class
# Returns accuracy, estimated classes, and probabilities for any single segment.
# segment.distances : should be a single matrix from the array of distances.
# classes : classification identifies for set being analyzed.
# k.grid : tuning parameteres to analyze.  Either neighbor size (integers)
#          or bandwidth constants (positive real values)
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above)

segment.class<-function(segment.distances, classes, k.grid=2, class.method="wknn", ker=kern.tri)
{
  n.objs<-dim(segment.distances)[1]
  closest.temp<-matrix(ncol=n.objs,nrow=n.objs-1)
  classify.temp<-matrix(ncol=n.objs,nrow=n.objs-1)
  dist.ord<-matrix(ncol=n.objs, nrow=n.objs-1)
  kern.probs<-matrix(ncol=n.objs, nrow=n.objs-1)
  
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  met.temp<-segment.distances
  
  for(j in 1:n.objs) 
  {
    or.temp<-order(met.temp[,j])
    closest.temp[,j]<-or.temp[or.temp!=j]
    classify.temp[,j]<-classes[closest.temp[,j]]
    dist.ord[,j]<-met.temp[or.temp[or.temp!=j], j]
  }
  
  if(class.method=="wknn")
  {
    kern.probs<-apply(dist.ord, 2, ker, h=1, method=class.method)
    prob.classes<-array(dim=c(n.objs, length(class.levels), length(k.grid)), dimnames = list(seq(1:n.objs), class.levels, seq(1:length(k.grid))))
    for(q in 1:length(k.grid))
    {
      k=k.grid[q]
      for(i in 1:length(class.levels.index))
      {
        for(j in 1:n.objs) 
        {
          assign.temp<-which(classify.temp[(1:k),j]==class.levels[i])
          total.prob<-sum(kern.probs[(1:k), j])  
          prob.classes[j,class.levels.index[i],q]=sum(kern.probs[assign.temp,j])/total.prob
        }
      }
    }
    
    est.classes<-matrix(nrow=n.objs, ncol=length(k.grid))
    for (k in 1:length(k.grid))
    {
      for (j in 1:n.objs)
      {
        max.classes<-class.levels[as.numeric(which(prob.classes[j,,k]==max(prob.classes[j,,k])))]
        est.classes[j,k]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
      }
    }
    
    est.accuracy<-numeric()
    for (k in 1:length(k.grid))
    {
      est.accuracy[k]<-mean(est.classes[,k]==classes)
    }
  } 
  
  if(class.method=="kernel")
  {
    prob.classes<-array(dim=c(n.objs, length(class.levels), length(k.grid)), dimnames = list(seq(1:n.objs), class.levels, seq(1:length(k.grid))))
    for(q in 1:length(k.grid))
    {
      kern.probs<-apply(dist.ord, 2, ker, h=k.grid[q], method=class.method)
      for(i in 1:length(class.levels.index))
      {
        for(j in 1:n.objs) 
        {
          assign.temp<-which(classify.temp[,j]==class.levels[i])
          total.prob<-sum(kern.probs[, j])
          ifelse(total.prob==0, 
                 prob.classes[j,class.levels.index[i],q]<-1/length(class.levels.index),
                 prob.classes[j,class.levels.index[i],q]<-sum(kern.probs[assign.temp,j])/total.prob
          )
        }
      }
    }
    
    est.classes<-matrix(nrow=n.objs, ncol=length(k.grid))
    for (k in 1:length(k.grid))
    {
      for (j in 1:n.objs)
      {
        max.classes<-class.levels[as.numeric(which(prob.classes[j,,k]==max(prob.classes[j,,k])))]
        est.classes[j,k]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
      }
    }
    
    est.accuracy<-numeric()
    for (k in 1:length(k.grid))
    {
      est.accuracy[k]<-mean(est.classes[,k]==classes)
    }
  }
  
  return(list(accuracy.est=est.accuracy, classes.est=est.classes, prob.array=prob.classes))
  
}

### grid.class
# Performs the standard grid search over tuning parameters and segments.
# No sequential stepping.  To be used if a predetermined grid of parameters is desired rather than using
# sequential stepping of segment sizes.
# distance.array : array of distances calculated from fdata object.  Distances for each segmented FDO.
# classes : classification identifies for set being analyzed.
# segments.grid : segment sizes to analyzes.  Should be vector of integers.
# k.grid : tuning parameteres to analyze.  Either neighbor size (integers)
#          or bandwidth constants (positive real values)
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above)
# step.method : stewpwise procedure.  Either forward.ensemble or all.segs.ensemble
# seg.weight : Set to TRUE to weight combination of segments by individual segment LOOCV accuracy. 
# thresh : stepwise accuracy improvement threshold.  Does not affect all.segs.ensemble.  
# do.par : Set to TRUE to run calculations in parallel.
# max.cores : number of cores to use during parallel calculations.
# output : Toggle for displaying information about calculation while calculation is running

grid.class<-function(distance.array, classes, segments.grid=1, k.grid=5, class.method="wknn", ker=kern.tri, step.method=forward.ensemble, seg.weight=TRUE, thresh=0.0001, do.par=FALSE, max.cores=2, output=FALSE)
{
  #if(output) cat("Evaluating", segments.grid, "Segment(s) \n")
  if(sum(segments.grid)!=dim(distance.array)[3]) stop("Segments and Distance Array do not match!")
  n.objs<-dim(distance.array)[1]
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  seg.start<-cumsum(segments.grid)-segments.grid
  
  if(output) cat("Evaluating Segment Probabilities. \n")
  
  s=1
  segment.array<-array(dim=c(n.objs, length(class.levels), length(k.grid), segments.grid[s]))
  accuracy.mat<-matrix(nrow=length(k.grid), ncol=segments.grid[s])
  est.classes<-matrix(nrow=n.objs, ncol=length(k.grid))
  
  if(do.par && segments.grid[s]!=1)
  {
    use.cores<-min(max.cores, segments.grid[s])
    cl.max<-makeCluster(use.cores)
    registerDoParallel(cl.max)
    kernel.output<-foreach(q=1:segments.grid[s], .packages = "fda.usc", .export="segment.class", .combine=append) %dopar%
    {
      return(list(segment.class(distance.array[,,seg.start[s]+q], classes, k.grid, class.method, ker)))
    }
    stopCluster(cl.max)
  }
  
  if(!do.par || segments.grid[s]==1)
  {
    kernel.output<-foreach(q=1:segments.grid[s], .packages = "fda.usc", .export="segment.class", .combine=append) %do%
    {
      return(list(segment.class(distance.array[,,seg.start[s]+q], classes, k.grid, class.method, ker)))
    }
  }
  
  for(q in 1:segments.grid[s])
  {
    segment.array[,,,q]<-kernel.output[[q]]$prob.array
    accuracy.mat[,q]<-kernel.output[[q]]$accuracy.est
  }
  
  if(output) cat("Finished Probabilities. \n")
  
  grid.results<-step.method(segment.array, accuracy.mat, classes, seg.weight=seg.weight, thresh=thresh, do.par=do.par, cores=max.cores)$ens.accuracies
  
  if(output) cat("Finished Ensembling. \n")
  
  c.names<-NULL
  for(j in 1:length(segments.grid)) c.names<-c(c.names, paste(segments.grid[j], "segs"))
  r.names<-NULL
  for(j in 1:length(k.grid)) r.names<-c(r.names, paste("k=", round(k.grid[j], 6), sep=""))
  grid.results<-as.matrix(grid.results)
  colnames(grid.results)<-c.names
  rownames(grid.results)<-r.names
  
  return(grid.results)
}

### seg.grid.class
# Evaluates a sequentially increasing segment size, evaluating if the top models (as determined by LOOCV)
# improve in accuracy as segment size increases. The top top.models.eval models must not change 
# for seg.sep segments.  Same functionality as grid.class but with improved selection of segment sizes.
# Can drastically improve computation times if segments to be analyzed is unknown. 
# Has additionally been improved to run directly from the fdata.object of interest.
# Indent and deriv must be supplied if fdata.object is to be manipulated within the calculation,
# or fdata.object can be supplied as derivative of interest and set deriv=0.
# fdata.object : set to be analyzed.  Default setting is original data curves, that can be differentiated
# if desired.
# classes : classification identifiers for set being analyzed.
# top.models : number of models to be considered when evaluating if segmentation has been optimized.
# seg.sep : distance between total segments evaluated and segment size of top models.  Setting to 0 will
#           stop calculation once min.segments has been analyzed.
# min.segments : initial segment sizes to be analyzed.  Should be a minimum of 2, although will work from 1.
# max.segments : stopping parameter to ensure sequential increasing of segment size will stop.
# indent : numerical noise indentation. 
# deriv : which derivative of fdata.object to be analyzed.    
# k.grid : grid of tuning parameters   
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above)
# step.method : stewpwise procedure.  Either forward.ensemble or all.segs.ensemble
# seg.weight : Set to TRUE to weight combination of segments by individual segment LOOCV accuracy
# thresh : stepwise accuracy improvement threshold.  Does not affect all.segs.ensemble.
# density : integration grid mesh size to control accuracy
# do.par : Set to TRUE to run calculations in parallel.
# max.cores : number of cores to use during parallel calculations.
# output : Toggle for displaying information about calculation while calculation is running
# write.out : Toggle to write data to file as calculation is performed.
# write.name : File name to which data is sent if write.out = T.

seq.grid.class<-function(fdata.object, classes, top.models.eval=15, seg.sep=1, min.segments=5, max.segments=30, 
                         indent=0, deriv=0, k.grid=c(1:10), class.method="wknn", ker=kern.tri, 
                         step.method=forward.ensemble, seg.weight=FALSE, thresh=0.0001, density=NULL,
                         do.par=FALSE, max.cores=2, output=FALSE, write.out=FALSE, write.name=NULL)
{
  if(is.null(write.name)) write.name="output"
  
  fd.full<-fdata.object
  argvals<-fdata.object$argvals
  if(deriv!=0) fd.full<-fdata.deriv(fd.full,deriv)
  n.objs<-length(classes)
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  full.grid<-NULL
  for(seg.size in 1:min.segments)
  {
    if(output) cat("Calculating Distance for", seg.size, "segment(s). \n")
    dist.temp<-calc.distance.array(fd.full, indent, density, total.segments=seg.size, do.par, max.cores)
    grid.seg<-grid.class(dist.temp, classes, segments.grid=seg.size, k.grid=k.grid, class.method=class.method, ker=ker, 
                         step.method=step.method, seg.weight=seg.weight, thresh=thresh, do.par=do.par, max.cores=max.cores, output=output)
    full.grid<-cbind(full.grid, grid.seg)  
  }
  
  if(write.out) write.csv(full.grid, paste(write.name, "_Grid.csv", sep=""))
  
  largest.seg.anal<-min.segments
  largest.seg.in.top<-max(as.numeric(substring(sort.models(full.grid, top.models.eval, do.print=FALSE)[1:top.models.eval,2], 1, nchar(sort.models(full.grid, top.models.eval, do.print=FALSE)[1:top.models.eval,2])-4)))
  
  while((largest.seg.in.top+seg.sep)>largest.seg.anal && min.segments<max.segments)
  {
    min.segments<-min.segments+1
    if(output) cat("Segment size increased to ", min.segments, ".\n")
    if(output) cat("Calculating Distance for", min.segments, "segment(s). \n")
    dist.temp<-calc.distance.array(fd.full, indent, density, total.segments=min.segments, do.par, max.cores)
    grid.seg<-grid.class(dist.temp, classes, segments.grid=min.segments, k.grid=k.grid, class.method=class.method, ker=ker, 
                         step.method=step.method, seg.weight=seg.weight, thresh=thresh, do.par=do.par, max.cores=max.cores, output=output)
    full.grid<-cbind(full.grid, grid.seg)
    if(write.out) write.csv(full.grid, paste(write.name, "_Grid.csv", sep=""))
    largest.seg.anal<-min.segments
    largest.seg.in.top<-max(as.numeric(substring(sort.models(full.grid, top.models.eval, do.print=FALSE)[1:top.models.eval,2], 1, nchar(sort.models(full.grid, top.models.eval, do.print=FALSE)[1:top.models.eval,2])-4)))
  }
  
  if(output) cat("Finished Sequential Grid Analysis at a total of", min.segments, "Segments. \n")
  return(full.grid)
}

### validation.probs
# Function that takes a truncated distance matrix (including only the testing data ) and a k.grid
# and returns a probability array for the test data.  Used for Cross validation.
# dist.temp : temporary truncated matrix of distances (distances from training to test objects)
# taining.classes : classification identifiers for training set
# k.eval : tuning parmeter
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above)

validation.probs<-function(dist.temp, training.classes, k.eval, class.method="wknn", ker=kern.tri)
{
  test.objs<-ncol(dist.temp)
  train.objs<-nrow(dist.temp)
  closest.temp<-matrix(ncol=test.objs, nrow=train.objs)
  classify.temp<-matrix(ncol=test.objs, nrow=train.objs)
  dist.ord<-matrix(ncol=test.objs, nrow=train.objs)
  kern.probs<-matrix(ncol=test.objs, nrow=train.objs)
  
  class.levels<-as.numeric(levels(factor(training.classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  for(j in 1:test.objs)
  {
    or.temp<-order(dist.temp[,j])
    closest.temp[,j]<-or.temp
    classify.temp[,j]<-training.classes[closest.temp[,j]]
    dist.ord[,j]<-dist.temp[or.temp,j]
  }
  
  if(class.method=="wknn")
  {
    kern.probs<-apply(dist.ord, 2, ker, h=1, method=class.method)
    prob.classes<-matrix(nrow=test.objs, ncol=length(class.levels))
    for(i in 1:length(class.levels.index))
    {
      for(j in 1:test.objs) 
      {
        assign.temp<-which(classify.temp[(1:k.eval),j]==class.levels[i])
        total.prob<-sum(kern.probs[(1:k.eval), j])  
        prob.classes[j,class.levels.index[i]]=sum(kern.probs[assign.temp,j])/total.prob
      }
    }
  }
  
  if(class.method=="kernel")
  {
    prob.classes<-matrix(nrow=test.objs, ncol=length(class.levels))
    kern.probs<-apply(dist.ord, 2, ker, h=k.eval, method=class.method)
    for(i in 1:length(class.levels.index))
    {
      for(j in 1:test.objs) 
      {
        assign.temp<-which(classify.temp[,j]==class.levels[i])
        total.prob<-sum(kern.probs[, j])
        ifelse(total.prob==0, 
               prob.classes[j,class.levels.index[i]]<-1/length(class.levels.index),
               prob.classes[j,class.levels.index[i]]<-sum(kern.probs[assign.temp,j])/total.prob
        )
      }
    }
  }
  
  return(prob.classes)
}

### grid.model.cv
# Runs cross validation of the top models.analyzed models from the results of a sequential grid search.
# fdata.object : set to be analyzed.  Default setting is original data curves.
# classes : classification identifies for set being analyzed.
# grid.results : output from running a grid search
# models.analyzed : number of models to validated
# folds : number of folds
# trials : number of times cross validation should be performed.
# folds.list : a list of folds identifiers
# seg.sep : distance between total segments evaluated and segment size of top models.  Setting to 0 will
#           stop calculation once min.segments has been analyzed.
# min.segments : initial segments sizes to be analyzed.  Should be a minimum of 2, although will work from 1.
# max.segments : stopping parameter to ensure sequential increasing of segment size will stop.
# indent : numerical noise indentation.
# deriv : which derivative of fdata.object to be analyzed.   
# k.grid : grid of tuning parameters   
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above)
# step.method : stewpwise procedure.  Either forward.ensemble or all.segs.ensemble
# seg.weight : Set to TRUE to weight combination of segments by individual segment LOOCV accuracy
# thresh : stepwise accuracy improvement threshold.  Does not affect all.segs.ensemble.
# density : integration grid mesh size to control accuracy 
# do.par : Set to TRUE to run calculations in parallel.
# max.cores : number of cores to use during parallel calculations.
# large.set : A toggle for numerical stability if analyzing a dataset with extreme number of obersvations
# output : Toggle for displaying information about calculation while calculation is running

grid.model.cv<-function(fdata.object, classes, grid.results, models.analyzed=10, folds=10, trials=1, 
                        folds.list=NULL, indent=0, deriv=0, density=NULL, class.method="wknn", ker=kern.tri,
                        step.method=forward.ensemble, seg.weight=FALSE, do.par=FALSE, max.cores=2, 
                        large.set=FALSE, output=TRUE)
{
  if(is.null(folds.list)) folds.list<-fold.creator(classes.temp, folds, trials)
  if(deriv!=0) fdata.object<-fdata.deriv(fdata.object, deriv)
  
  new.k<-as.numeric(substring(sort.models(grid.results, models.analyzed, do.print=FALSE)[1:models.analyzed,1], 3))
  new.segs<-as.numeric(substring(sort.models(grid.results, models.analyzed, do.print=FALSE)[1:models.analyzed,2], 1, nchar(sort.models(grid.results, models.analyzed, do.print=FALSE)[1:models.analyzed,2])-4))
  seg.grid<-as.numeric(names(table(new.segs)))
  
  n.objs<-length(classes)
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  train.all.acc.mat<-matrix(nrow=trials*folds, ncol=models.analyzed)
  test.all.acc.mat<-matrix(nrow=trials*folds, ncol=models.analyzed)
  models.out<-list()
  
  for(seg.chosen in 1:length(seg.grid))
  {
    total.segments=seg.grid[seg.chosen]
    models.to.analyze<-which(new.segs==total.segments)
    model.k<-new.k[models.to.analyze]
    
    if(output) cat("Segment", seg.chosen, "of ", length(seg.grid), "\n")
    
    distance.array<-calc.distance.array(fdata.object, indent, density, total.segments=total.segments, do.par, max.cores)
    
    prob.array.temp<-array(dim=c(n.objs, length(class.levels.index), length(model.k), total.segments))
    acc.mat<-matrix(ncol=length(model.k), nrow=total.segments)
    
    for(j in 1:total.segments)
    {
      temp.1<-segment.class(distance.array[,,j], classes, model.k, class.method, ker)
      prob.array.temp[,,,j]<-temp.1$prob.array
      acc.mat[j,]<-temp.1$accuracy.est
    }
    
    model.temp<-step.method(prob.array.temp, t(acc.mat), classes, seg.weight=seg.weight)
    model.segs.used<-model.temp$ens.segments
    models.out[[seg.chosen]]<-model.segs.used
    
    if(do.par && !large.set)
    {
      cl.red<-makeCluster(max.cores) 
      registerDoParallel(cl.red)
      cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %dopar%
      {
        est.accuracy.train<-matrix(nrow=folds, ncol=length(model.k))
        est.accuracy.test<-matrix(nrow=folds, ncol=length(model.k))
        
        for(k.fold in 1:folds)
        {
          
          test.folds<-which(folds.list[[j.trial]]==k.fold)
          training.classes<-classes[-test.folds]
          test.classes<-classes[test.folds]
          test.objs=length(test.classes)
          training.objs<-length(training.classes)
          
          # dist.temp.test<-distance.array[-test.folds,test.folds,]
          # dist.temp.train<-distance.array[-test.folds,-test.folds,]
          
          for(model in 1:length(model.k))
          {
            len.model.segs<-length(model.segs.used[[model]])
            
            ### Determine Training Set Accuracy and Segment LOO Accuracy for each Training Set ###    
            prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
            acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
            
            for(j in 1:len.model.segs)
            {
              temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                    training.classes, model.k[model], class.method, ker)
              prob.array.train[,,1,j]<-temp.1$prob.array
              acc.mat.train[j,]<-temp.1$accuracy.est
            }
            
            ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
            for(q in class.levels.index)
            {
              ifelse(seg.weight,
                     ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                     ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
              )
              if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
            }
            
            est.classes.train<-numeric()
            for(q in 1:training.objs)
            {
              max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
              est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
            }
            
            est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)

            prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
            
            for(q in 1:len.model.segs)
              prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model], class.method, ker)                 
            
            ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
            for(q in class.levels.index)
            {
              ifelse(seg.weight,
                     ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                     ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
              )
              if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
            }
            
            est.classes.test<-numeric()
            for(q in 1:test.objs)
            {
              max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
              est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
            }
            
            est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
            
          }
        }
        
        return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
        
      }
      stopCluster(cl.red)
    }
    
    
    if(!do.par || large.set)
    {
      cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %do%
      {
        est.accuracy.train<-matrix(nrow=folds, ncol=length(model.k))
        est.accuracy.test<-matrix(nrow=folds, ncol=length(model.k))
        
        for(k.fold in 1:folds)
        {
          
          test.folds<-which(folds.list[[j.trial]]==k.fold)
          training.classes<-classes[-test.folds]
          test.classes<-classes[test.folds]
          test.objs=length(test.classes)
          training.objs<-length(training.classes)
          
          for(model in 1:length(model.k))
          {
            len.model.segs<-length(model.segs.used[[model]])
            prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
            acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
            
            for(j in 1:len.model.segs)
            {
              temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                    training.classes, model.k[model], class.method, ker)
              prob.array.train[,,1,j]<-temp.1$prob.array
              acc.mat.train[j,]<-temp.1$accuracy.est
            }
            
            ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
            for(q in class.levels.index)
            {
              ifelse(seg.weight,
                     ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                     ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
              )
              if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
            }
            
            est.classes.train<-numeric()
            for(q in 1:training.objs)
            {
              max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
              est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
            }
            
            est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)

            prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
            
            for(q in 1:len.model.segs)
              prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model], class.method, ker)                 
            
            ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
            for(q in class.levels.index)
            {
              ifelse(seg.weight,
                     ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                     ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
              )
              if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
            }
            
            est.classes.test<-numeric()
            for(q in 1:test.objs)
            {
              max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
              est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
            }
            
            est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
            
          }
        }
        
        return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
        
      }
    }
    
    test.accs.temp<-NULL
    train.accs.temp<-NULL
    for(j in 1:trials)
    {
      test.accs.temp<-rbind(test.accs.temp, cv.out[[j]]$test.accuracy)
      train.accs.temp<-rbind(train.accs.temp, cv.out[[j]]$train.accuracy)
    }
    
    for(j in 1:length(models.to.analyze))
    {
      train.all.acc.mat[,models.to.analyze[j]]<-train.accs.temp[,j]
      test.all.acc.mat[,models.to.analyze[j]]<-test.accs.temp[,j]
    }
  }
  
  new.comb<-cbind(new.k, new.segs)
  x<-NULL
  for(j in 1:length(models.out)) x<-c(x, lengths(models.out[[j]]))
  segments.used.grid<-matrix(NA, ncol=models.analyzed, nrow=max(x))
  
  for(j in 1:length(seg.grid))
  {
    col.temp<-which(new.comb[,2]==seg.grid[[j]])
    n.temp<-length(col.temp)
    for(k in 1:n.temp)
    {
      segs.used<-models.out[[j]][[k]]
      len.temp<-length(segs.used)
      segments.used.grid[(1:len.temp),col.temp]<-segs.used
    }
  }
  
  seg.sum<-rbind(t(new.comb), segments.used.grid)
  
  return(list(test.accuracies=test.all.acc.mat, training.accuracies=train.all.acc.mat, segment.summary=seg.sum))
}

### dist.to.prob
# Takes a distance array and produces probability array for given method, kernel, and turning parameter
# See segment.class 
# distance.array : array of distances 
# total.segments : total number of segments (should match array dimension) 
# k.size : a single tuning parameter (not a vector!) 
# classes : classification identifies for set being analyzed. Only used for dimensionality. 

dist.to.prob<-function(distance.array, total.segments, k.size, classes)
{
  n.objs<-length(classes)
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  temp.array<-array(dim=c(n.objs,n.objs,total.segments))
  temp.array[,,1:total.segments]<-distance.array
  
  prob.array.temp<-array(dim=c(n.objs, length(class.levels.index), 1, total.segments))
  acc.mat<-matrix(ncol=1, nrow=total.segments)
  
  for(j in 1:total.segments)
  {
    temp.1<-segment.class(temp.array[,,j], classes.temp, k.size, class.method, ker)
    prob.array.temp[,,,j]<-temp.1$prob.array
    acc.mat[j,]<-temp.1$accuracy.est
  }
  
  return(list(probability.array=prob.array.temp, accuracies=acc.mat))
}

### bestsub.ensemble
# Takes a probability array and evaluates all combinations reporting the top LOOCV models for each combination size
# step.array : array of probabilities determined from distances, and dependent on which method and kernel chosen.
#              Calculated frequently in other functions.
# segment.asccuracies : accuracy of each segment included in step.array.  Used for weighted combinations.
# classes : classification identifies for set being analyzed.
# seg.weight : Set to TRUE if combining segment probabilities should be weighted by individual segment accuracy
# do.par : Set to TRUE to run calculations in parallel.
# max.cores : number of cores to use during parallel calculations.
# best.sub.max : maximum combination size to analyze (computational stability)

bestsub.ensemble<-function(step.array, segment.accuracies, classes, seg.weight=FALSE, do.par=FALSE, 
                           max.cores=2, best.sub.max=10)
{
  n.objs<-dim(step.array)[1]
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  knn.n<-dim(step.array)[3]
  segment.n<-dim(step.array)[4]
  seg.seq<-seq(1,segment.n,1)
  
  ens.acc<-numeric()
  ens.segs.used<-list()
  probs.ens<-array(dim=c(n.objs, length(class.levels.index), knn.n))
  classes.ens<-matrix(nrow=n.objs, ncol=knn.n)
  
  k=1
  
  max.vals<-which(segment.accuracies[k,]==max(segment.accuracies[k,]))
  max.seg<-ifelse(length(max.vals)==1, max.vals, sample(max.vals,1))
  
  segs.best<-max.seg
  acc.best<-segment.accuracies[k,max.seg]
  probs.best<-step.array[,,k,max.seg]
  
  est.classes.best<-numeric()
  for(q in 1:n.objs)
  {
    max.classes<-class.levels[as.numeric(which(probs.best[q,]==max(probs.best[q,])))]
    est.classes.best[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1)) #if tie sample from mean distances from the test objects. Mean distance over all segments?
  }
  
  best.sub.start<-list(acc=acc.best, segs=segs.best, probs=probs.best, class=est.classes.best)
  
  ifelse(segment.n>best.sub.max, use.segs.tot<-best.sub.max, use.segs.tot<-segment.n)
  
  if(do.par && segment.n!=1)
  {
    cl.max<-makeCluster(max.cores)
    registerDoParallel(cl.max)
    best.sub.list<-foreach(m=2:use.segs.tot, .combine=append) %dopar%
    {
      combs<-combn(segment.n,m)
      best.sub.temp<-list()
      best.acc=0
      
      for(j in 1:ncol(combs))
      {
        segs.update<-combs[,j]
        probs.temp<-matrix(nrow=n.objs, ncol=length(class.levels.index))
        for(q in class.levels.index)
        {
          ifelse(seg.weight,
                 probs.temp[,q]<-rowSums(t(segment.accuracies[k,segs.update]*t(as.matrix(step.array[,q,k,segs.update])))/sum(segment.accuracies[k,segs.update])),
                 probs.temp[,q]<-rowSums(as.matrix(step.array[,q,k,segs.update]))/length(segs.update)
          )
        }
        
        est.classes.temp<-numeric()
        for(q in 1:n.objs)
        {
          max.classes<-class.levels[as.numeric(which(probs.temp[q,]==max(probs.temp[q,])))]
          est.classes.temp[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
        }
        
        acc.temp<-mean(est.classes.temp==classes)
        
        if(acc.temp==best.acc) #acc.temp==acc.best
        {
          best.sub.temp<-append(best.sub.temp, list(list(acc=acc.temp, segs=segs.update, probs=probs.temp, class=est.classes.temp)))
        }
        
        if(acc.temp>best.acc)
        {
          best.acc=acc.temp
          best.sub.temp<-list(list(acc=acc.temp, segs=segs.update, probs=probs.temp, class=est.classes.temp))
        }
        
        
      }
      
      return(list(best.sub.temp))
    }
    
    stopCluster(cl.max)  
  }
  
  if(!do.par && segment.n!=1)
  {
    best.sub.list<-foreach(m=2:use.segs.tot, .combine=append) %do%
    {
      combs<-combn(segment.n,m)
      best.sub.temp<-list()
      best.acc=0
      
      for(j in 1:ncol(combs))
      {
        segs.update<-combs[,j]
        probs.temp<-matrix(nrow=n.objs, ncol=length(class.levels.index))
        for(q in class.levels.index)
        {
          ifelse(seg.weight,
                 probs.temp[,q]<-rowSums(t(segment.accuracies[k,segs.update]*t(as.matrix(step.array[,q,k,segs.update])))/sum(segment.accuracies[k,segs.update])),
                 probs.temp[,q]<-rowSums(as.matrix(step.array[,q,k,segs.update]))/length(segs.update)
          )
        }
        
        est.classes.temp<-numeric()
        for(q in 1:n.objs)
        {
          max.classes<-class.levels[as.numeric(which(probs.temp[q,]==max(probs.temp[q,])))]
          est.classes.temp[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
        }
        
        acc.temp<-mean(est.classes.temp==classes)
        
        if(acc.temp==best.acc) #acc.temp==acc.best
        {
          best.sub.temp<-append(best.sub.temp, list(list(acc=acc.temp, segs=segs.update, probs=probs.temp, class=est.classes.temp)))
        }
        
        if(acc.temp>best.acc)
        {
          best.acc=acc.temp
          best.sub.temp<-list(list(acc=acc.temp, segs=segs.update, probs=probs.temp, class=est.classes.temp))
        }
        
        
      }
      
      return(list(best.sub.temp))
    }
  }
  
  acc.tab<-numeric()
  
  n.top.models=1
  if(segment.n!=1)
  {
    for(k in 1:(use.segs.tot-1))
    {
      for(j in 1:length(best.sub.list[[k]]))
      {
        n.top.models=n.top.models+1
      }  
    }  
  }

  segs.tab<-matrix(NA, ncol=n.top.models, nrow=use.segs.tot)
  acc.tab[1]<-best.sub.start$acc
  segs.tab[1,1]<-best.sub.start$segs
  count=1
  
  if(segment.n!=1)
  {
    for(k in 1:(use.segs.tot-1))
    {
      for(j in 1:length(best.sub.list[[k]]))
      {
        count=count+1
        acc.tab[count]<-best.sub.list[[k]][[j]]$acc
        segs.tab[1:(k+1),count]<-best.sub.list[[k]][[j]]$segs
      }  
    }  
  }

  return(list(accuracies=acc.tab, segments.used=segs.tab))
  
}

###  bss.single.model.cv
#  Evaluates all best segment selection specifications for a given model (chosen segment size and neighbor size)

bss.single.model.cv<-function(fdata.object, classes, k.size=5, total.segments=1, 
                              indent=0, deriv=0, best.sub.max=10, density=NULL, folds=10, 
                              trials=1, folds.list=NULL, smooth=FALSE, class.method="wknn", 
                              ker=kern.tri, seg.weight=FALSE, do.par=FALSE, max.cores=2, large.set=FALSE, output=FALSE)
{
  if(output) time.start<-Sys.time()
  n.objs<-length(classes)
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  if(is.null(folds.list)) folds.list<-fold.creator(classes.temp, folds, trials)
  if(deriv!=0) fdata.object<-fdata.deriv(fdata.object, deriv)
  
  distance.array<-calc.distance.array(fdata.object, indent, density, total.segments=total.segments, do.par, max.cores)
  
  if(output) 
  {
    time.dist<-Sys.time()
    cat("Finished Distances.", time.dist-time.start, "\n")
  }
  
  temp.prob<-dist.to.prob(distance.array, total.segments, k.size, classes.temp)
  
  temp<-bestsub.ensemble(step.array=temp.prob$probability.array, segment.accuracies=t(temp.prob$accuracies), 
                         classes=classes.temp, seg.weight=seg.weight, do.par=do.par, 
                         max.cores=max.cores, best.sub.max=best.sub.max)
  
  if(output) 
  {
    time.bss<-Sys.time()
    cat("Finished Best Segment Selection", time.bss-time.dist, "\n")
  }
  
  segs.list<-list()
  for(j in 1:ncol(temp[[2]]))
  {
    segs.list[[j]]<-as.numeric(na.exclude(temp[[2]][,j]))
  }
  
  model.segs.used<-segs.list
  models.to.analyze<-1:length(model.segs.used)
  model.k<-rep(k.size, length(model.segs.used))
  train.all.acc.mat<-matrix(nrow=trials*folds, ncol=length(model.segs.used))
  test.all.acc.mat<-matrix(nrow=trials*folds, ncol=length(model.segs.used))
  
  if(do.par && !large.set)
  {
    cl.red<-makeCluster(max.cores)
    registerDoParallel(cl.red)
    cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %dopar%
    {
      est.accuracy.train<-matrix(nrow=folds, ncol=length(model.k))
      est.accuracy.test<-matrix(nrow=folds, ncol=length(model.k))
      
      for(k.fold in 1:folds)
      {
        
        test.folds<-which(folds.list[[j.trial]]==k.fold)
        training.classes<-classes[-test.folds]
        test.classes<-classes[test.folds]
        test.objs=length(test.classes)
        training.objs<-length(training.classes)
        
        # dist.temp.test<-distance.array[-test.folds,test.folds,]
        # dist.temp.train<-distance.array[-test.folds,-test.folds,]
        
        for(model in 1:length(model.k))
        {
          len.model.segs<-length(model.segs.used[[model]])
          
          ### Determine Training Set Accuracy and Segment LOO Accuracy for each Training Set ###    
          prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
          acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
          
          for(j in 1:len.model.segs)
          {
            temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                  training.classes, model.k[model], class.method, ker)
            prob.array.train[,,1,j]<-temp.1$prob.array
            acc.mat.train[j,]<-temp.1$accuracy.est
          }
          
          ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
          }
          
          est.classes.train<-numeric()
          for(q in 1:training.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
            est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)
          
          ### Determine Test Set Accuracy ###
          
          prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
          
          for(q in 1:len.model.segs)
            prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model], class.method, ker)                 
          
          ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
          }
          
          est.classes.test<-numeric()
          for(q in 1:test.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
            est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
          
        }
      }
      
      return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
      
      
      est.accuracy.test
      
    }
    stopCluster(cl.red)
  }
  
  if(!do.par || large.set)
  {
    cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %do%
    {
      est.accuracy.train<-matrix(nrow=folds, ncol=length(model.k))
      est.accuracy.test<-matrix(nrow=folds, ncol=length(model.k))
      
      for(k.fold in 1:folds)
      {
        
        test.folds<-which(folds.list[[j.trial]]==k.fold)
        training.classes<-classes[-test.folds]
        test.classes<-classes[test.folds]
        test.objs=length(test.classes)
        training.objs<-length(training.classes)
        
        # dist.temp.test<-distance.array[-test.folds,test.folds,]
        # dist.temp.train<-distance.array[-test.folds,-test.folds,]
        
        for(model in 1:length(model.k))
        {
          len.model.segs<-length(model.segs.used[[model]])
          
          ### Determine Training Set Accuracy and Segment LOO Accuracy for each Training Set ###    
          prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
          acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
          
          for(j in 1:len.model.segs)
          {
            temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                  training.classes, model.k[model], class.method, ker)
            prob.array.train[,,1,j]<-temp.1$prob.array
            acc.mat.train[j,]<-temp.1$accuracy.est
          }
          
          ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
          }
          
          est.classes.train<-numeric()
          for(q in 1:training.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
            est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)
          
          ### Determine Test Set Accuracy ###
          
          prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
          
          for(q in 1:len.model.segs)
            prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model], class.method, ker)                 
          
          ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
          }
          
          est.classes.test<-numeric()
          for(q in 1:test.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
            est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
          
        }
      }
      
      return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
      
    }
  }
  
  test.accs.temp<-NULL
  train.accs.temp<-NULL
  for(j in 1:trials)
  {
    test.accs.temp<-rbind(test.accs.temp, cv.out[[j]]$test.accuracy)
    train.accs.temp<-rbind(train.accs.temp, cv.out[[j]]$train.accuracy)
  }
  
  for(j in 1:length(models.to.analyze))
  {
    train.all.acc.mat[,models.to.analyze[j]]<-train.accs.temp[,j]
    test.all.acc.mat[,models.to.analyze[j]]<-test.accs.temp[,j]
  }
  
  return(list(test.accuracies=test.all.acc.mat, training.accuracies=train.all.acc.mat, bss=temp))
}

### bss.ens.model.cv
# Cross validates the top LOOCV accuracy segment combinations for all combination sizes using
# classifiers produced from a single curve or combination of multiple derivative orders.
# If a single curve is to be analyzed, k.size, seg.size, indent, deriv should be single elements.
# If multiple curves are to be combined, k.size, seg.size, indent, deriv should be vectors.
# best.sub.max can be used to control the size of the segment combinations evaluated.
# fdata.object : original data functional object.
# classes : classification identifies for set being analyzed.
# k.sizes : chosen model tuning parameters for each curve to be analyzed.
# seg.size : chosen model segment sizes for each curve to be analyzed.
# indent : numerical derivative noise indents for each derivative order.
# deriv : derivative orders.  i.e c(0,1,2) would be original curve with first and second derivs.
# best.sub.max : maximum combination size.
# thresh : only if forward segment selection it to be used.  Improvement in LOOCV to add segment to ensemble.
# density : integration grid mesh size to control accuracy  
# folds : number of folds.
# trials : number of times cross validation should be performed.
# folds.list : a list of folds identifiers
# class.method : nonparametric classification method.  Either "wknn" or "kernel".
# ker : kernel to be used (see Kernels above).
# seg.weight : Set to TRUE if combining segment probabilities should be weighted by individual segment accuracy.
# use.forward : Use FSS in place of BSS when evaluating final model.  Reduces computational burden but
#               does not guarantee the top segment combination will be found.
# do.par : Set to TRUE to run calculations in parallel.
# max.cores : number of cores to use during parallel calculations.
# large.set : A toggle for numerical stability if analyzing a dataset with extreme number of obersvations
# output : Toggle for displaying information about calculation while calculation is running

bss.ens.model.cv<-function(fdata.object, classes, k.sizes=c(5,5), seg.sizes=c(1,2),
                           indent=c(0,0), deriv=c(0,1), best.sub.max=10, thresh=1e-4, density=NULL, folds=10,
                           trials=1, folds.list=NULL, class.method="wknn", ker=kern.tri, seg.weight=FALSE,
                           use.forward=FALSE, do.par=FALSE, max.cores=2, large.set=FALSE, output=FALSE)
{
  if(output) time.start<-Sys.time()
  n.objs<-length(classes)
  class.levels<-as.numeric(levels(factor(classes)))
  class.levels.index<-seq(1:length(class.levels))
  
  if(is.null(folds.list)) folds.list<-fold.creator(classes.temp, folds, trials)
  total.segments<-sum(seg.sizes)
  
  distance.array<-array(dim=c(n.objs, n.objs, total.segments))
  temp.prob.array<-array(dim=c(n.objs, length(class.levels.index), 1, total.segments))
  temp.prob.accur<-NULL
  
  csum<-cumsum(seg.sizes)
  csum.seg<-(cumsum(seg.sizes)-seg.sizes)+1
  
  for(j in 1:length(seg.sizes))
  {
    ifelse(deriv[j]==0, fdata.object.temp<-fdata.object, fdata.object.temp<-fdata.deriv(fdata.object, deriv[j]))
    distance.array[,,csum.seg[j]:csum[j]]<-calc.distance.array(fdata.object.temp, indent[j], density, total.segments=seg.sizes[j], do.par, max.cores)
    temp.prob<-dist.to.prob(distance.array[,,csum.seg[j]:csum[j]], seg.sizes[j], k.sizes[j], classes.temp)
    temp.prob.array[,,,csum.seg[j]:csum[j]]<-temp.prob$probability.array
    temp.prob.accur<-c(temp.prob.accur, temp.prob$accuracies)
  }
  
  if(output)
  {
    time.dist<-Sys.time()
    cat("Finished Distances.", time.dist-time.start, "\n")
  }
  
  if(!use.forward)
  {
    temp<-bestsub.ensemble(step.array=temp.prob.array, segment.accuracies=t(as.matrix(temp.prob.accur)), 
                           classes=classes.temp, seg.weight=seg.weight, do.par=do.par, 
                           max.cores=max.cores, best.sub.max=best.sub.max)
    
    segs.list<-list()
    for(j in 1:ncol(temp[[2]]))
    {
      segs.list[[j]]<-as.numeric(na.exclude(temp[[2]][,j]))
    }
    
    model.segs.used<-segs.list
    models.to.analyze<-1:length(model.segs.used)
  }
  
  if(use.forward)
  {
    temp.fss<-forward.ensemble(step.array=temp.prob.array, segment.accuracies=t(as.matrix(temp.prob.accur)), 
                               classes=classes.temp, seg.weight=seg.weight, thresh=thresh, do.par=do.par, 
                               cores=max.cores)
    
    temp<-list(segments.used=as.matrix(temp.fss$ens.segments[[1]], ncol=1))
    
    model.segs.used<-list(temp.fss$ens.segments[[1]])
    models.to.analyze<-1
  }
  
  if(output)
  {
    time.bss<-Sys.time()
    cat("Finished BSS.", time.bss-time.dist, "\n")
  }
  
  model.k<-NULL
  deriv.k<-NULL
  for(j in 1:length(seg.sizes)) 
  {
    deriv.k<-c(deriv.k, rep(deriv[j], seg.sizes[j]))
    model.k<-c(model.k, rep(k.sizes[j], seg.sizes[j]))
  }
  
  train.all.acc.mat<-matrix(nrow=trials*folds, ncol=length(model.segs.used))
  test.all.acc.mat<-matrix(nrow=trials*folds, ncol=length(model.segs.used))
  
  if(do.par && !large.set)
  {
    cl.red<-makeCluster(max.cores)
    registerDoParallel(cl.red)
    cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %dopar%
    {
      est.accuracy.train<-matrix(nrow=folds, ncol=length(model.segs.used))
      est.accuracy.test<-matrix(nrow=folds, ncol=length(model.segs.used))
      
      for(k.fold in 1:folds)
      {
        test.folds<-which(folds.list[[j.trial]]==k.fold)
        training.classes<-classes[-test.folds]
        test.classes<-classes[test.folds]
        test.objs=length(test.classes)
        training.objs<-length(training.classes)
        
        for(model in 1:length(model.segs.used))
        {
          len.model.segs<-length(model.segs.used[[model]])
          prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
          acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
          
          for(j in 1:len.model.segs)
          {
            temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                  training.classes, model.k[model.segs.used[[model]][j]], class.method, ker)
            prob.array.train[,,1,j]<-temp.1$prob.array
            acc.mat.train[j,]<-temp.1$accuracy.est
          }
          
          ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
          }
          
          est.classes.train<-numeric()
          for(q in 1:training.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
            est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)

          prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
          
          for(q in 1:len.model.segs)
            prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model.segs.used[[model]][q]], class.method, ker)                 
          
          ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
          }
          
          est.classes.test<-numeric()
          for(q in 1:test.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
            est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
          
        }
      }
      
      return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
      
    }
    stopCluster(cl.red)
  }
  
  if(!do.par || large.set)
  {
    cv.out<-foreach(j.trial=1:trials, .packages=c("fda.usc"),.export=c("segment.class", "validation.probs", "kern.tri", "kern.norm", "kern.unif"),.combine=append) %do%
    {
      est.accuracy.train<-matrix(nrow=folds, ncol=length(model.segs.used))
      est.accuracy.test<-matrix(nrow=folds, ncol=length(model.segs.used))
      
      for(k.fold in 1:folds)
      {
        test.folds<-which(folds.list[[j.trial]]==k.fold)
        training.classes<-classes[-test.folds]
        test.classes<-classes[test.folds]
        test.objs=length(test.classes)
        training.objs<-length(training.classes)
        
        for(model in 1:length(model.segs.used))
        {
          len.model.segs<-length(model.segs.used[[model]])
          prob.array.train<-array(dim=c(training.objs, length(class.levels.index), 1, len.model.segs))
          acc.mat.train<-matrix(ncol=1, nrow=len.model.segs)
          
          for(j in 1:len.model.segs)
          {
            temp.1<-segment.class(distance.array[-test.folds,-test.folds, model.segs.used[[model]][j]], 
                                  training.classes, model.k[model.segs.used[[model]][j]], class.method, ker)
            prob.array.train[,,1,j]<-temp.1$prob.array
            acc.mat.train[j,]<-temp.1$accuracy.est
          }
          
          ens.probs.train<-matrix(nrow=training.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.train[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.train[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.train[,q]<-rowSums(as.matrix(prob.array.train[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.train[,q])) probs.temp[,q]<-rep(0, train.objs)
          }
          
          est.classes.train<-numeric()
          for(q in 1:training.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.train[q,]==max(ens.probs.train[q,])))]
            est.classes.train[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.train[k.fold, model]<-mean(est.classes.train==training.classes)

          prob.array.test<-array(dim=c(test.objs, length(class.levels.index), 1, len.model.segs))
          
          for(q in 1:len.model.segs)
            prob.array.test[,,1,q]<-validation.probs(distance.array[-test.folds,test.folds,model.segs.used[[model]][q]], training.classes, k.eval=model.k[model.segs.used[[model]][q]], class.method, ker)                 
          
          ens.probs.test<-matrix(nrow=test.objs, ncol=length(class.levels.index))
          for(q in class.levels.index)
          {
            ifelse(seg.weight,
                   ens.probs.test[,q]<-rowSums(t(acc.mat.train[1:len.model.segs,]*t(as.matrix(prob.array.test[,q,,1:len.model.segs]))))/sum(acc.mat.train[1:len.model.segs,]),
                   ens.probs.test[,q]<-rowSums(as.matrix(prob.array.test[,q,,1:len.model.segs]))/len.model.segs
            )
            if(anyNA(ens.probs.test[,q])) probs.temp[,q]<-rep(0, test.objs)
          }
          
          est.classes.test<-numeric()
          for(q in 1:test.objs)
          {
            max.classes<-class.levels[as.numeric(which(ens.probs.test[q,]==max(ens.probs.test[q,])))]
            est.classes.test[q]<-ifelse(length(max.classes)==1, max.classes, sample(max.classes,1))
          }
          
          est.accuracy.test[k.fold, model]<-mean(est.classes.test==test.classes)
          
        }
      }
      
      return(list(list(test.accuracy=est.accuracy.test, train.accuracy=est.accuracy.train)))
      
    }
  }
  
  test.accs.temp<-NULL
  train.accs.temp<-NULL
  for(j in 1:trials)
  {
    test.accs.temp<-rbind(test.accs.temp, cv.out[[j]]$test.accuracy)
    train.accs.temp<-rbind(train.accs.temp, cv.out[[j]]$train.accuracy)
  }
  
  for(j in 1:length(models.to.analyze))
  {
    train.all.acc.mat[,models.to.analyze[j]]<-train.accs.temp[,j]
    test.all.acc.mat[,models.to.analyze[j]]<-test.accs.temp[,j]
  }
  
  return(list(test.accuracies=test.all.acc.mat, training.accuracies=train.all.acc.mat, bss=temp))
}