###
# Ensemble of Segmented Functional Nonparametric Classifiers
# Robert Buscaglia (rbuscagl@asu.edu)
# Yiannis Kamarianakis (yiannis76@asu.edu)
# June 20, 2017
#
# This is supplementary code supporting the submitted manuscript "Ensemble of Segmented
# Functional Nonparametric Classifiers".  Included in the repository along
# with this file are the following.
#
#
# ESFuNC_Functions.R : Contains all required functions and packages.
# tecator.txt : Text file containing the 215 Tecator samples. First 215 columns are data
#               points, column 216 is the class (0=lowfat, 1=highfat). data points were collected 
#               at wavelength values from 850 to 1048 nm with spacing of 2 nm.
# lupustherm.txt : Text file containing the 589 Lupus Thermogram samples.  First 451 columns
#                  are data points, column 452 is the class (0=nonlupus, 1=lupus).  Data points were
#                  collected at temperatures of 45 to 90 every 0.1 celsius.
#
#
# This template is meant to aid the user in reproducing the Tecator dataset results.  The 
# Lupus results can also be analyzed, but the user will be required to make modifications to the
# file. Comments or questions should be emailed to the authors above.
#
# The template will export several files.  For the code below to work correctly, some of the
# exported files will be loaded from disk.  This is done to both allow the user to see
# all of the results of running a particular section, while also allowing the sections to be
# run at different times.  Although the Tecator dataset can be analyzed in a short time frame,
# for larget data sets the user may wish to run the algorithm in steps.  Exporting all of the
# results ensures that no information is lost.
#
# I would encourage exploring the code as there are a lot of details that go into the analysis.
# Be aware though that trying different classification methods ("wknn" or "kernel"), different kernels,
# or different datasets will require the user to manipulate the code below.
#
# If one is only interested in obtaining the Tecator results, the file can be sourced after 
# changing the working directory in Step 1.  For reference, using 4 of the 8 available cores
# of an Intel Core i7-6700 takes ~8 minutes.


######################################
### Step 1 - Set working directory ###
######################################
# This should be the location on computer containing all of the files from the Github repository.
# The algorithm generates a significant amount of output.  This is also the location to which all
# files will be exported.  Defaults to C:\\

setwd("C:\\")


##############################################
### Step 2 - Source functions and packages ###
##############################################
# Sources the Functions and Packages required to run the algorithm

source("ESFuNC_Functions.R")


#############################
### Step 3 - Load dataset ###
#############################
# Template is for Tecator data set analysis.
# Required from this step is that data.temp be the original data points collected at argvals.
# Also defines class.temp, a vector of classification identifies. (0 = lowfat, 1=highfat for Tecator)
# The Tecator data points were collected at wavelength values from 850 to 1048 nm with spacing of 2 nm.

tecator<-read.table("tecator.txt")
data.temp<-as.matrix(tecator[,1:100])
classes.temp<-as.numeric(tecator[,101])
argvals=seq(850, 1048, 2)
rm(tecator)


#######################################################
### Step 4 - Create fdata object from original data ###
#######################################################
# By default creates B-spline functional representations of the original data.
# The impact of smoothing, basis representation, and data reduction are being investigated.
# This step requires that the functional data object be created for the original data, but
# alterations to how the object is created should be accomidated.

##### smooth=FALSE

fdata.temp<-fdata(data.temp, argvals)


#################################################
### Step 5 - Data Set and Analysis Parameters ###
#################################################
# This sets the parameters of the analysis.  Each parameter has a brief description.
#
# Data Set Information

data.set="tecator_default" #name holder of dataset being analyzed
result.name = "wknn_tri" #name holder of nonparametric method being used

indent.set=c(0,0,0) #indentation of dataset to remove numerical  noise.  None required for Tecator.
                    #LUPUS NOTE: indent.set=c(0,3,6)
deriv.set=c(0,1,2)  #derivatives to analyze.  deriv=0 is equivalent to analyzing original dataset.

# Grid Parameters

min.segments=6 #minimum number of segments that should be analyzed during grid search
max.segments=30 #maximim number of segments.  Used to ensure algorithm stops.
top.models.eval=15 #number of top LOOCV specifications to evaluate.
seg.sep=2 #seperation between segment sizes necessary to halt the sequential grid search
k.grid<-1:100 #neighborhood size, when using wknn must be integers
# k.grid<-exp(seq(log(1e-4),log(10),(log(10)-log(1e-4))/(250-1))) #a potential grid if using class.method="kernel"
class.method="wknn" #classification method, can be either "wknn" or "kernel.  Default "wknn"
ker=kern.tri #kernel to be used when generating probabilities.  See KERNELS in ESFuNC_Functions.R
step.method=forward.ensemble #Stepwise method to be used.  Only forward.ensemble or all.segs.ensemble currently available.
seg.weight=TRUE #weighting segment combinations by LOOCV accuracy individual segments.
thresh=0.0001 #stepwise improvement in LOOCV necassary to include segment in model.
density=NULL #integration density.  NULL (Recommended) will use original grid mesh density for all integration steps. 
do.par=TRUE #use parallelization when possible.  True (Recommended)
max.cores=4 #number of cores to use during parallelization.  Default=4. Have used up to 22 on server without issue.
output=FALSE #option to monitor computations.  Not fully functional but can be useful during long computations.
write.out=FALSE #option to export files as computations is running.  May not be fully functional.
write.name.temp=NULL #name to use if exporting files during computation.

# Cross Validation Parameters

models.analyzed=30 #number of specifications from grid search to evaluate
folds=10 #number of folds to use
trials=20 #number of times K-fold cross validation is repeated.
best.sub.max=10

## Folds
# This creates a default list of fold identifiers so that all cross validation steps are performed using the
# same folds.  Strongly Recommended, but any list of folds can be used.

set.seed(1)
folds.list<-fold.creator(classes.temp, folds, trials) #fold identifiers for multiple runs (trials) of K-fold CV


######################################################################################
### Step 6 - Evaluate Sequential Grid Search and Cross Validate Top Specifications ###
######################################################################################
# This section performs the sequential grid search for reach derivative order, exporting the resulting
# grid to a .csv file.  The grid search results are then used to evaluate the top performing specfications.
# This includes a segmentation size and tuning parameter.  From cross validation, a model is selected
# with preference given to more parsimonious models.

###
# Original Data
###

write.name.0=paste(data.set,"_",result.name,"_0", sep="") #Name Holder

##run grid search

curve.0.grid<-seq.grid.class(fdata.object=fdata.temp, classes=classes.temp, top.models.eval=top.models.eval, 
                          seg.sep=seg.sep, min.segments=min.segments, max.segments=max.segments, 
                          indent=indent.set[1], deriv=deriv.set[1], k.grid=k.grid, class.method=class.method, ker=ker, 
                          step.method=step.method, seg.weight=seg.weight, thresh=thresh, density=density,
                          do.par=do.par, max.cores=max.cores, output=output, write.out=write.out, write.name=write.name.temp)

write.csv(curve.0.grid, paste(write.name.0, "_Grid.csv", sep="")) #export Grid
sort.models(curve.0.grid, models.analyzed) #View top specifications

##cross validate top specifications

curve.0.grid.cv<-grid.model.cv(fdata.object=fdata.temp, classes=classes.temp, grid.results=curve.0.grid, 
                               models.analyzed=models.analyzed, folds=folds, trials=trials, folds.list=folds.list,
                               indent=indent.set[1], deriv=deriv.set[1], density=density, class.method=class.method,
                               ker=ker, step.method=step.method, seg.weight=seg.weight, do.par=do.par, max.cores=max.cores,
                               large.set=FALSE, output=output)

##export grid results and create image

write.csv(curve.0.grid.cv$test.accuracies, paste(write.name.0, "_test_CV.csv", sep=""))
write.csv(curve.0.grid.cv$training.accuracies, paste(write.name.0, "_train_CV.csv", sep=""))
write.csv(curve.0.grid.cv$segment.summary, paste(write.name.0, "_segment-summary_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))

###
# First Derivative
###

write.name.1=paste(data.set,"_",result.name,"_1", sep="") #Name Holder

##run grid search

curve.1.grid<-seq.grid.class(fdata.object=fdata.temp, classes=classes.temp, top.models.eval=top.models.eval, 
                             seg.sep=seg.sep, min.segments=min.segments, max.segments=max.segments, 
                             indent=indent.set[2], deriv=deriv.set[2], k.grid=k.grid, class.method=class.method, ker=ker, 
                             step.method=step.method, seg.weight=seg.weight, thresh=thresh, density=density,
                             do.par=do.par, max.cores=max.cores, output=output, write.out=write.out, write.name=write.name.temp)

write.csv(curve.1.grid, paste(write.name.1, "_Grid.csv", sep=""))
sort.models(curve.1.grid, models.analyzed)

##cross validate top specifications

curve.1.grid.cv<-grid.model.cv(fdata.object=fdata.temp, classes=classes.temp, grid.results=curve.1.grid, 
                               models.analyzed=models.analyzed, folds=folds, trials=trials, folds.list=folds.list,
                               indent=indent.set[2], deriv=deriv.set[2], density=density, class.method=class.method,
                               ker=ker, step.method=step.method, seg.weight=seg.weight, do.par=do.par, max.cores=max.cores,
                               large.set=FALSE, output=output)

##export grid results and create image

write.csv(curve.1.grid.cv$test.accuracies, paste(write.name.1, "_test_CV.csv", sep=""))
write.csv(curve.1.grid.cv$training.accuracies, paste(write.name.1, "_train_CV.csv", sep=""))
write.csv(curve.1.grid.cv$segment.summary, paste(write.name.1, "_segment-summary_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))

###
# Second Derivative
###

write.name.2=paste(data.set,"_",result.name,"_2", sep="") #Name Holder

##run grid sesarch

curve.2.grid<-seq.grid.class(fdata.object=fdata.temp, classes=classes.temp, top.models.eval=top.models.eval, 
                             seg.sep=seg.sep, min.segments=min.segments, max.segments=max.segments, 
                             indent=indent.set[3], deriv=deriv.set[3], k.grid=k.grid, class.method=class.method, ker=ker, 
                             step.method=step.method, seg.weight=seg.weight, thresh=thresh, density=density,
                             do.par=do.par, max.cores=max.cores, output=output, write.out=write.out, write.name=write.name.temp)



write.csv(curve.2.grid, paste(write.name.2, "_Grid.csv", sep=""))
sort.models(curve.2.grid, models.analyzed)

##cross validate top specifications

curve.2.grid.cv<-grid.model.cv(fdata.object=fdata.temp, classes=classes.temp, grid.results=curve.2.grid, 
                               models.analyzed=models.analyzed, folds=folds, trials=trials, folds.list=folds.list,
                               indent=indent.set[3], deriv=deriv.set[3], density=density, class.method=class.method,
                               ker=ker, step.method=step.method, seg.weight=FALSE, do.par=do.par, max.cores=max.cores,
                               large.set=FALSE, output=output)

##export grid results and create image

write.csv(curve.2.grid.cv$test.accuracies, paste(write.name.2, "_test_CV.csv", sep=""))
write.csv(curve.2.grid.cv$training.accuracies, paste(write.name.2, "_train_CV.csv", sep=""))
write.csv(curve.2.grid.cv$segment.summary, paste(write.name.2, "_segment-summary_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))


###
#                                              PAUSE
# NOTE : As long as datafiles have been exported, calculations can be paused here without issues.  All
# required data will be loaded from disk when necassary.  If you choose to pause here, always ensure
# that Steps 1 and 2 have been run before continuing and that the image of previous steps has been saved.
###


#######################################################
### Step 7 - Choose Model for Each Derivative Order ###
#######################################################
# The results of cross validation are used to choose a top performing model for each derivative order.
# The CV results are read from disk.  The models are compared using mean test set accuracy.  All specifications
# that are within one standard error of the maximum mean test set accuracy are considered equivalent.
# If model models return equivalent mean test set accuracy, the model with the smallest segmentation size
# and largest tuning parameter is chosen for parsimony.
#
# This sections also makes for a good moment to evaluate how the algorithm is working.  It is possible
# to make boxplots of the cross-validation results, more detailed summaries, or investiage the results
# further if desired.  Some code for boxplots is shown below, and the create.cv.sum gives median and MAD
# statistics and can be viewed easily.
#
# Model selection has not been automated, and the model must be chosen by the user.  The correct models
# for the Tecator dataset have been entered.  If no changes have been made to the file (except working directory)
# then the user should not have to make any decisions below.
# 

###
# Read CV Results from Disk
###

load(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))
cv.test.0<-read.csv(paste(write.name.0, "_test_CV.csv", sep=""), row.names=1)
cv.test.1<-read.csv(paste(write.name.1, "_test_CV.csv", sep=""), row.names=1)
cv.test.2<-read.csv(paste(write.name.2, "_test_CV.csv", sep=""), row.names=1)

###
# Create summary statistics for each model
###

cv.test.0.sum<-create.cv.sum(cv.test.0)
cv.test.1.sum<-create.cv.sum(cv.test.1)
cv.test.2.sum<-create.cv.sum(cv.test.2)

###
# Original Data - Choose a Model
###

which(cv.test.0.sum[4,]>=max(cv.test.0.sum[4,])-cv.test.0.sum[5,]/sqrt(trials*folds))
#sort.models(curve.0.grid, models.analyzed)
#boxplot(cv.test.0)

model.0.chosen=1
k.0<-as.numeric(curve.0.grid.cv$segment.summary[1,model.0.chosen])
s.0<-as.numeric(curve.0.grid.cv$segment.summary[2,model.0.chosen])

###
# First Derivative - Choose a Model
###

which(cv.test.1.sum[4,]>=max(cv.test.1.sum[4,])-cv.test.1.sum[5,]/sqrt(trials*folds))
#sort.models(curve.1.grid, models.analyzed)
#boxplot(cv.test.1)

model.1.chosen=1
k.1<-as.numeric(curve.1.grid.cv$segment.summary[1,model.1.chosen])
s.1<-as.numeric(curve.1.grid.cv$segment.summary[2,model.1.chosen])

###
# Second Derivative - Choose a Model
###

which(cv.test.2.sum[4,]>=max(cv.test.2.sum[4,])-cv.test.2.sum[5,]/sqrt(trials*folds))
sort.models(curve.2.grid, models.analyzed)
#boxplot(cv.test.2)

model.2.chosen=1
k.2<-as.numeric(curve.2.grid.cv$segment.summary[1,model.2.chosen])
s.2<-as.numeric(curve.2.grid.cv$segment.summary[2,model.2.chosen])

save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))


#######################################################
### Step 8 - Evaluate Best Segment Selection Models ###
#######################################################
# This section takes the models chosen above and evaluates best segment selection results using the given
# tuning parameter and segmentation size for each derivative order. he individual curve models (original,
# first, and second derivatives) are evaluated.  The final model is produced by evaluating the combination
# of models.  Any combination can be evaluated (i.e. original + first derivative, original + second derivative).
# The only combination coded below is for the resuls of combining all three curves (original + first derivative +
# second derivative).  The user can easily manipuate the given functions if other combinations are of interest.
#
# This function (which can be found in the ESFuNC_Functions.R file) works by first determining the top LOOCV
# best segment combination model for each segment size up to the given best.sub.max.  Each of the top models
# is then evaluated by k-fold cross validation.  The model which returns the highest mean test set accuracy
# is considered the final model.  Users are encouraged to evaluate the output in detail, as many interesting values
# are exported.  This includes how the segment combination changes as combination size is increased.  It is possible
# to find jumps that would  not be possible to make when using forward segment selection (i.e. the top model for
# including one more segment may remove two of the old segments, which would not be possible using FSS).  The
# results of BSS and cross validation of BSS models is very interesting.  Future work using machine learning 
# techniques may aide these steps greatly.
#
# For the Tecator dataset, the results of combining all three curves is done quickly.  Caution should be taken
# when evaluating resulting models that give large segmentation sizes.  When the combined total segments reaches
# 20 or more segments, computations will slow drastically due to BSS growing exponentially.  At 30 or more segments
# the computations may have significant trouble completing unless using a significant number of cores.  The user
# should be warned not to try to combine 30 or more segments using BSS.
#
# A setting has been added that will allow for final models to be determined from Forward Segment selection if desired.
# The user may enter the option : use.forward = TRUE
#
# Once the bss.ens.model.cv has been run, the resulting solutions can be found in the .csv files.

###
# Original Curves BSS and CV
###

bss.0.cv<-bss.ens.model.cv(fdata.object=fdata.temp, classes=classes.temp, k.sizes=c(k.0), seg.sizes=c(s.0),
                           indent=indent.set[c(1)], deriv=deriv.set[c(1)], density=density, folds=folds, trials=trials,
                           folds.list=folds.list, class.method=class.method, ker=ker, seg.weight=seg.weight,
                           do.par=do.par, max.cores=max.cores, large.set=FALSE, output=output, best.sub.max=10)

write.csv(bss.0.cv$test.accuracies, paste(write.name.0, "_test_BSS_CV.csv", sep=""))
write.csv(bss.0.cv$training.accuracies, paste(write.name.0, "_train_BSS_CV.csv", sep=""))
write.csv(rbind(bss.0.cv$bss$accuracies, bss.0.cv$bss$segments.used), paste(write.name.0, "_segment-summary_BSS_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))

###
# First Derivative BSS and CV
###

bss.1.cv<-bss.ens.model.cv(fdata.object=fdata.temp, classes=classes.temp, k.sizes=c(k.1), seg.sizes=c(s.1),
                           indent=indent.set[c(2)], deriv=deriv.set[c(2)], density=density, folds=folds, trials=trials,
                           folds.list=folds.list, class.method=class.method, ker=ker, seg.weight=seg.weight,
                           do.par=do.par, max.cores=max.cores, large.set=FALSE, output=output, best.sub.max=10)

write.csv(bss.1.cv$test.accuracies, paste(write.name.1, "_test_BSS_CV.csv", sep=""))
write.csv(bss.1.cv$training.accuracies, paste(write.name.1, "_train_BSS_CV.csv", sep=""))
write.csv(rbind(bss.1.cv$bss$accuracies, bss.1.cv$bss$segments.used), paste(write.name.1, "_segment-summary_BSS_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))

###
# Second Derivative BSS and CV
###

bss.2.cv<-bss.ens.model.cv(fdata.object=fdata.temp, classes=classes.temp, k.sizes=c(k.2), seg.sizes=c(s.2),
                           indent=indent.set[c(3)], deriv=deriv.set[c(3)], density=density, folds=folds, trials=trials,
                           folds.list=folds.list, class.method=class.method, ker=ker, seg.weight=seg.weight,
                           do.par=do.par, max.cores=max.cores, large.set=FALSE, output=output, best.sub.max=10)

write.csv(bss.2.cv$test.accuracies, paste(write.name.2, "_test_BSS_CV.csv", sep=""))
write.csv(bss.2.cv$training.accuracies, paste(write.name.2, "_train_BSS_CV.csv", sep=""))
write.csv(rbind(bss.2.cv$bss$accuracies, bss.2.cv$bss$segments.used), paste(write.name.2, "_segment-summary_BSS_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))


###
# Final Model - All Curves Combined BSS and CV
###

write.name.0.1.2=paste(data.set,"_",result.name,"_0_1_2", sep="")

bss.0.1.2.cv<-bss.ens.model.cv(fdata.object=fdata.temp, classes=classes.temp, k.sizes=c(k.0,k.1,k.2), seg.sizes=c(s.0,s.1,s.2),
                           indent=indent.set[c(1,2,3)], deriv=deriv.set[c(1,2,3)], density=density, folds=folds, trials=trials,
                           folds.list=folds.list, class.method=class.method, ker=ker, seg.weight=seg.weight,
                           do.par=do.par, max.cores=max.cores, large.set=FALSE, output=output, best.sub.max=10)

write.csv(bss.0.1.2.cv$test.accuracies, paste(write.name.0.1.2, "_test_BSS_CV.csv", sep=""))
write.csv(bss.0.1.2.cv$training.accuracies, paste(write.name.0.1.2, "_train_BSS_CV.csv", sep=""))
write.csv(rbind(bss.0.1.2.cv$bss$accuracies, bss.0.1.2.cv$bss$segments.used), paste(write.name.0.1.2, "_segment-summary_BSS_CV.csv", sep=""))
save.image(paste(data.set,"_",result.name,"_IMAGE.RData", sep=""))


##################################################
### Step 9 - Visualization of Top ESFuNC Model ###
##################################################
# This section recreates the figure from the manuscript, and should be flexible to any other results or data sets
# that are analyzed.  Image is exported to disk as .tiff of publication quality.
#
# To reproduce the figure exactly, cross validation of functional weighted-knn without any segmetation is performed.

load(paste(data.set,"_",result.name,"_IMAGE.RData", sep="")) #Load Image.

### Perform functional knn cross validation (no segmentation)
# 

new.mat<-as.matrix(curve.0.grid[,1])
colnames(new.mat)<-"1 segs"

single.segment.cv<-grid.model.cv(fdata.object=fdata.temp, classes=classes.temp, grid.results=new.mat, 
                                 models.analyzed=10, folds=folds, trials=trials, folds.list=folds.list,
                                 indent=indent.set[1], deriv=deriv.set[1], density=density, class.method=class.method,
                                 ker=ker, step.method=step.method, seg.weight=seg.weight, do.par=do.par, max.cores=max.cores,
                                 large.set=FALSE, output=output)


single.seg.model<-which.max(colMeans(single.segment.cv$test.accuracies))

#boxplot(single.segment.cv$test.accuracies)
write.csv(single.segment.cv$test.accuracies, paste(data.set,"_",result.name,"_Single_Segment_CV.csv", sep=""))

### Setup Functional Objects for Plotting
# 

fdata.temp.0<-fdata.temp
fdata.temp.1<-fdata.deriv(fdata.temp, deriv.set[2])
fdata.temp.2<-fdata.deriv(fdata.temp, deriv.set[3])

### Set up segment boundaries
# This is how partitioning is done in create.segments, but we need everything explicit here for viewing.

min.0<-argvals[1+indent.set[1]]
max.0<-argvals[length(argvals)-indent.set[1]]
step.0<-(max.0-min.0)/(s.0)
splits.0<-numeric()
for(j in 1:(s.0+1)) splits.0[j]<-min.0+(j-1)*step.0

min.1<-argvals[1+indent.set[2]]
max.1<-argvals[length(argvals)-indent.set[2]]
step.1<-(max.1-min.1)/(s.1)
splits.1<-numeric()
for(j in 1:(s.1+1)) splits.1[j]<-min.1+(j-1)*step.1

min.2<-argvals[1+indent.set[3]]
max.2<-argvals[length(argvals)-indent.set[3]]
step.2<-(max.2-min.2)/(s.2)
splits.2<-numeric()
for(j in 1:(s.2+1)) splits.2[j]<-min.2+(j-1)*step.2

### Setup Boxplots for Comparison
# Boxplots compare how well ESFuNC performs on each individual derivative order compared to 
# ensembling segments from multiple derivative orders.

test.cv.0.bss<-read.csv(paste(data.set,"_", result.name, "_0_test_BSS_cv.csv", sep=""), row.names=1)
test.0.bss.sum<-create.cv.sum(test.cv.0.bss)
test.0.bss.chosen<-which.max(test.0.bss.sum[4,])

test.cv.1.bss<-read.csv(paste(data.set,"_", result.name, "_1_test_BSS_cv.csv", sep=""), row.names=1)
test.1.bss.sum<-create.cv.sum(test.cv.1.bss)
test.1.bss.chosen<-which.max(test.1.bss.sum[4,])

test.cv.2.bss<-read.csv(paste(data.set,"_", result.name, "_2_test_BSS_cv.csv", sep=""), row.names=1)
test.2.bss.sum<-create.cv.sum(test.cv.2.bss)
test.2.bss.chosen<-which.max(test.2.bss.sum[4,])

### Setup Final BSS Model from Multiple Derivative Orders
# This also obtains the segments that are used in the final model, which is used when plotting.

seg.sum.0.1.2.bss<-read.csv(paste(data.set,"_", result.name, "_0_1_2_segment-summary_BSS_CV.csv", sep=""))
seg.sum.0.1.2.bss<-round(seg.sum.0.1.2.bss[,-1],4)
test.cv.0.1.2.bss<-read.csv(paste(data.set,"_", result.name, "_0_1_2_test_BSS_cv.csv", sep=""), row.names=1)

test.0.1.2.bss.sum<-create.cv.sum(test.cv.0.1.2.bss)
test.0.1.2.bss.chosen<-which.max(test.0.1.2.bss.sum[4,])
segs.0.1.2<-c(1:s.0, 1:s.1, 1:s.2)
segs.inc.0.1.2<-na.omit(seg.sum.0.1.2.bss[-1,test.0.1.2.bss.chosen])
segs.exc.0.1.2<-segs.0.1.2[-segs.inc.0.1.2]

### Final Model Boxplots
#

final.model.box<-cbind(single.segment.cv$test.accuracies[,single.seg.model], test.cv.0.bss[,test.0.bss.chosen], test.cv.1.bss[,test.1.bss.chosen], test.cv.2.bss[,test.2.bss.chosen], test.cv.0.1.2.bss[,test.0.1.2.bss.chosen])

### Produce the Figure
# Everything has been formatting for the Tecator Figure. Using other data will require reformatting the code.

labels<-c(expression(paste("F"^"(0)", "-FKNN", sep="")), expression(paste("F"^"(0)", "-ESFNC", sep="")), expression(paste("F"^"(1)", "-ESFNC", sep="")), expression(paste("F"^"(2)", "-ESFNC", sep="")), expression(paste("F"^"(0,1,2)", "-ESFNC", sep="")))

tiff("Default_Figure.tiff", width=20, height=20, units="cm", res=600)

m<-rbind(c(1,2), c(3,4))
layout(m)

par(mar=c(4.0,4.5,0.5,0.5))
plot(fdata.temp.0, col="grey50", lty=1,
     main="", xlab="", ylab="", cex.axis=1.1, cex.lab=1.35, las=1, font.lab=2, yaxt='n')
axis(side=2, at=seq(-2,5,1), labels=T, las=1, cex.axis=1.1)
title(xlab="t", line=2.3, cex.lab=1.35, font.lab=2)
title(ylab="X(t)", line=2.5, cex.lab=1.35, font.lab=2)
abline(v=splits.0)
abline(h=0)
rect(-100,-100,splits.0[1],100, col="black")
rect(splits.0[length(splits.0)],-100,splits.0[length(splits.0)]+100,100, col="black")
for(j in setdiff(1:s.0,intersect(segs.inc.0.1.2,segs.0.1.2[1:(s.0)])))
  rect(splits.0[j],-100,splits.0[j+1],100, col="black", density=10)
mtext("(A)", side=2, line=2.5, cex=1.25, at=5.5, las=1)

par(mar=c(4.0,4.5,0.5,0.5))
plot(fdata.temp.1, col="grey50", lty=1, main="", xlab="", ylab="",
     cex.axis=1.1, cex.lab=1.35, las=1, font.lab=2, yaxt='n', ylim=c(-0.02,0.05))
axis(side=2, at=seq(-0.02, 0.05, 0.02), labels=T, las=1, cex.axis=1.1)
title(xlab="t", line=2.3, cex.lab=1.35, font.lab=2)
title(ylab="X'(t)", line=3, cex.lab=1.35, font.lab=2)
abline(v=splits.1)
abline(h=0)
rect(-100,-100,splits.1[1],10, col="black")
rect(splits.1[length(splits.1)],-10,splits.1[length(splits.1)]+100,10, col="black")
for(j in setdiff(1:s.1,intersect(segs.inc.0.1.2-s.0,segs.0.1.2[(s.0+1):(s.0+s.1)])))
  rect(splits.1[j],-10,splits.1[j+1],100, col="black", density=10)
mtext("(B)", side=2, line=3, cex=1.25, at=0.0505, las=1)

par(mar=c(3.5,4.5,1,0.5))
plot(fdata.temp.2, col="grey50", lty=1, main="", xlab="", ylab="",
     cex.axis=1.1, cex.lab=1.35, las=1, font.lab=2, ylim=c(-0.004, 0.004), yaxt='n')
axis(side=2, at=c(-0.004, -0.002, 0.002, 0.004), labels=T, las=1, cex.axis=0.9)
axis(side=2, at=c(0), labels=T, las=1, cex.axis=0.9)
title(xlab="t", line=2.3, cex.lab=1.35, font.lab=2)
title(ylab="X''(t)", line=2.5, cex.lab=1.35, font.lab=2)
abline(v=splits.2)
abline(h=0)
rect(-10,-10,splits.2[1],10, col="black")
rect(splits.2[length(splits.2)],-10,splits.2[length(splits.2)]+100,10, col="black")
for(j in setdiff(1:s.2,intersect(segs.inc.0.1.2-s.0-s.1,segs.0.1.2[(s.0+s.1+1):(s.0+s.1+s.2)])))
  rect(splits.2[j],-10,splits.2[j+1],100, col="black", density=10)
mtext("(C)", side=2, line=2.5, cex=1.25, at=0.005, las=1)

boxplot(final.model.box, las=2, ylab="", xlab="", xaxt='n', yaxt='n', cex.axis=1.1, cex.lab=1.25)
axis(side=1, at=seq(1,5,1), labels=FALSE, las=1, cex.axis=1.1)
text(seq(0.6,4.6,by=1), par("usr")[3]-0.025, labels=labels, srt=30, pos=1, xpd=TRUE, cex=0.9)
axis(side=2, at=seq(0,1,0.1), labels=T, las=1, cex.axis=1.1)
#title(xlab="Model", line=2.3, cex.lab=1.35, font.lab=2)
title(ylab="Accuracy", line=3.1, cex.lab=1.35, font.lab=2)
mtext("(D)", side=2, line=2.5, cex=1.25, at=1.04, las=1)

dev.off()
