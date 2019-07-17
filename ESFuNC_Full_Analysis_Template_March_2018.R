##########################################################
# ESFuNC Template Analysis
# Provides solutions using all three algorithm schemes.
# Dr. Robert Buscaglia
# February 2018
##########################################################

########################################################################
### Step 1 - Set working directory and Generalized Algorithm Settings 
########################################################################
# This should be the location on computer containing all of the files from the Github repository.
# The algorithm generates a significant amount of output.  This is also the location to which all
# files will be exported.  Defaults to C:\\

setwd("C:\\Users\\rober\\Dropbox\\Arizona State University\\Ioannis Research\\Fall 2017 Projects\\Prospectus\\R Code and Figures\\ESFuNC\\Results") #WD
source("ESFuNC_Functions_March_2018.R") #Updated_Functions

#####################
### Import dataset
#####################

cat("Setting Parameters. \n")

# lupuskernelnorm<-read.table("lupustherm.txt")
# data.temp<-as.matrix(lupuskernelnorm[,1:451])
# classes.temp<-as.numeric(lupuskernelnorm[,452])
# argvals=seq(45, 90, 0.1)

tecator<-read.table("tecator.txt")
data.temp<-as.matrix(tecator[,1:100])
classes.temp<-as.numeric(tecator[,101])
argvals=seq(850, 1048, 2)
rm(tecator)

output.name<-"TECATOR_WKNN_UNIF.RDATA"

###################################################################################
### Create basic fdata object : raw data discretization
### Give corresponding information about basis smoothing and roughness penalties
###################################################################################

fdata.temp<-fdata(data.temp, argvals=argvals)
nbasis.in<-100 ### Full Discretization
fpen.lambda.in<-0 ### No roughness

# gcv.basis<-61
# lambda.opt<-0.038
# out<-min.basis(fdata.temp, numbasis = 61, lambda=lambda.opt)
# plot(out$fdata.est)

################################################
### Derivatives to Analyze and Numerica Noise
################################################

deriv.analyze<-c(0,1,2)  ### Document Assumes use of original curves with first and second derivatives
indent.deriv<-c(0,3,6)

##################
### CV Settings
##################

folds=10
trials=20
set.seed(1)
folds.list=fold.creator(classes.temp, folds, trials)

#########################
### Parallel Computing
#########################

do.par=TRUE
max.cores=6

##########################
### Classifier Settings
##########################

class.method.in="wknn"
ker.in=kern.unif
k.grid.in<-2:100
density.in<-NULL

#################################
### Stepwise Ensemble Settings
#################################

step.method.in<-forward.ensemble
seg.weight.in<-TRUE
thresh.in<-0.0001

#############################
### Grid Search Parameters
#############################

top.models.eval.in<-15
seg.sep.in<-2
min.segments.in<-6
max.segments.in<-30

#############################
### Grid Validation Parameters
#############################

models.analyzed.in<-2*top.models.eval.in
best.sub.max.in<-25

###################################
### Hierarchical Search Settings
###################################

curve.seg.limit.single.in=20
curve.seg.limit.combined.in=30

######################################################################################################################################################
######################################################################################################################################################


#########################################################################
### Step 2 - Greedy Search !!!  Expensive! 
########################################################################



##############################
### Run Gridsearch on FDO 1
##############################

cat("grid.search.d0 \n")

grid.search.d0<-seq.grid.class(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, top.models.eval.in, seg.sep.in, min.segments.in, max.segments.in, 
                               indent.deriv[1], deriv.analyze[1], k.grid.in, class.method.in, ker.in, step.method.in, seg.weight.in, thresh.in, density.in,
                                do.par, max.cores, output=TRUE, write.out=FALSE, write.name=NULL)

save.image(output.name)

##############################
### Validate FDO 1 Top Models
##############################

cat("grid.search.d0.cv \n")

grid.search.d0.cv<-grid.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, grid.search.d0, models.analyzed.in, folds, trials, 
                              folds.list, indent.deriv[1], deriv.analyze[1], density.in, class.method.in, ker.in,
                              step.method.in, seg.weight.in, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

##############################
### Run Gridsearch on FDO 2
##############################

cat("grid.search.d1 \n")

grid.search.d1<-seq.grid.class(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, top.models.eval.in, seg.sep.in, min.segments.in, max.segments.in, 
                               indent.deriv[2], deriv.analyze[2], k.grid.in, class.method.in, ker.in, step.method.in, seg.weight.in, thresh.in, density.in,
                               do.par, max.cores, output=TRUE, write.out=FALSE, write.name=NULL)

save.image(output.name)

##############################
### Validate FDO 2 Top Models
##############################

cat("grid.search.d1.cv \n")

grid.search.d1.cv<-grid.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, grid.search.d0, models.analyzed.in, folds, trials, 
                                 folds.list, indent.deriv[2], deriv.analyze[2], density.in, class.method.in, ker.in,
                                 step.method.in, seg.weight.in, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

##############################
### Run Gridsearch on FDO 3
##############################

cat("grid.search.d2 \n")

grid.search.d2<-seq.grid.class(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, top.models.eval.in, seg.sep.in, min.segments.in, max.segments.in, 
                               indent.deriv[3], deriv.analyze[3], k.grid.in, class.method.in, ker.in, step.method.in, seg.weight.in, thresh.in, density.in,
                               do.par, max.cores, output=TRUE, write.out=FALSE, write.name=NULL)

save.image(output.name)

##############################
### Validate FDO 3 Top Models
##############################

cat("grid.search.d2.cv \n")

grid.search.d2.cv<-grid.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, grid.search.d0, models.analyzed.in, folds, trials, 
                                 folds.list, indent.deriv[3], deriv.analyze[3], density.in, class.method.in, ker.in,
                                 step.method.in, seg.weight.in, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

#########################################################################
### Set Top Performing Grid Search Models 
########################################################################

cat("Setting Top Models! \n")

d0.top.greedy.model<-which.max(colMeans(grid.search.d0.cv$test.accuracies))
d1.top.greedy.model<-which.max(colMeans(grid.search.d1.cv$test.accuracies))
d2.top.greedy.model<-which.max(colMeans(grid.search.d2.cv$test.accuracies))

k.d0<-as.numeric(grid.search.d0.cv$segment.summary[1,d0.top.greedy.model])
s.d0<-as.numeric(grid.search.d0.cv$segment.summary[2,d0.top.greedy.model])

k.d1<-as.numeric(grid.search.d1.cv$segment.summary[1,d1.top.greedy.model])
s.d1<-as.numeric(grid.search.d1.cv$segment.summary[2,d1.top.greedy.model])

k.d2<-as.numeric(grid.search.d2.cv$segment.summary[1,d2.top.greedy.model])
s.d2<-as.numeric(grid.search.d2.cv$segment.summary[2,d2.top.greedy.model])

save.image(output.name)


#####################################
### BSS Results of FDO 1 Top Model
#####################################

cat("bss.greedy.d0 \n")

ifelse(s.d0<25, use.forward.bss.d0<-FALSE, use.forward.bss.d0<-TRUE)

bss.greedy.d0<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=k.d0, seg.sizes=s.d0,
                                indent.deriv[1], deriv.analyze[1], best.sub.max.in, thresh.in, density.in, folds,
                                trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                use.forward.bss.d0, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

#####################################
### BSS Results of FDO 2 Top Model
#####################################

cat("bss.greedy.d1 \n")

ifelse(s.d1<25, use.forward.bss.d1<-FALSE, use.forward.bss.d1<-TRUE)

bss.greedy.d1<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=k.d1, seg.sizes=s.d1,
                                indent.deriv[2], deriv.analyze[2], best.sub.max.in, thresh.in, density.in, folds,
                                trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                use.forward.bss.d1, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

#####################################
### BSS Results of FDO 3 Top Model
#####################################

cat("bss.greedy.d2 \n")

ifelse(s.d2<25, use.forward.bss.d2<-FALSE, use.forward.bss.d2<-TRUE)

bss.greedy.d2<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=k.d2, seg.sizes=s.d2,
                                indent.deriv[3], deriv.analyze[3], best.sub.max.in, thresh.in, density.in, folds,
                                trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                use.forward.bss.d2, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

################################################
### Combined BSS Ensembles : All Combinations
################################################

##############################################
### BSS Results of FDO 1 + FDO 2 Top Models
##############################################

cat("bss.greedy.d0.d1 \n")

ifelse(sum(c(s.d0,s.d1))<30, use.forward.bss.d0.d1<-FALSE, use.forward.bss.d0.d1<-TRUE)

bss.greedy.d0.d1<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=c(k.d0,k.d1), seg.sizes=c(s.d0,s.d1),
                                indent.deriv[c(1,2)], deriv.analyze[c(1,2)], best.sub.max.in, thresh.in, density.in, folds,
                                trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                use.forward.bss.d0.d1, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

##############################################
### BSS Results of FDO 1 + FDO 3 Top Models
##############################################

cat("bss.greedy.d0.d2 \n")

ifelse(sum(c(s.d0,s.d2))<30, use.forward.bss.d0.d2<-FALSE, use.forward.bss.d0.d2<-TRUE)

bss.greedy.d0.d2<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=c(k.d0,k.d2), seg.sizes=c(s.d0,s.d2),
                                   indent.deriv[c(1,3)], deriv.analyze[c(1,3)], best.sub.max.in, thresh.in, density.in, folds,
                                   trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                   use.forward.bss.d0.d2, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

##############################################
### BSS Results of FDO 2 + FDO 3 Top Models
##############################################

cat("bss.greedy.d1.d2 \n")

ifelse(sum(c(s.d1,s.d2))<30, use.forward.bss.d1.d2<-FALSE, use.forward.bss.d1.d2<-TRUE)

bss.greedy.d1.d2<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=c(k.d1,k.d2), seg.sizes=c(s.d1,s.d2),
                                   indent.deriv[c(2,3)], deriv.analyze[c(2,3)], best.sub.max.in, thresh.in, density.in, folds,
                                   trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                   use.forward.bss.d1.d2, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

######################################################
### BSS Results of FDO 1 + FDO 2 + FDO 3 Top Models
######################################################

cat("bss.greedy.d0.d1.d2 \n")

ifelse(sum(c(s.d0,s.d1,s.d2))<30, use.forward.bss.d0.d1.d2<-FALSE, use.forward.bss.d0.d1.d2<-TRUE)

bss.greedy.d0.d1.d2<-bss.ens.model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, k.sizes=c(k.d0,k.d1,k.d2), seg.sizes=c(s.d0,s.d1,s.d2),
                                   indent.deriv[c(1,2,3)], deriv.analyze[c(1,2,3)], best.sub.max.in, thresh.in, density.in, folds,
                                   trials, folds.list, class.method.in, ker.in, seg.weight.in,
                                   use.forward.bss.d0.d1.d2, do.par, max.cores, large.set=FALSE, output=TRUE)

save.image(output.name)

######################################################################################################################################################
######################################################################################################################################################
### END GREEDY ANALYSIS - bss.greedy.* will contain all results.
######################################################################################################################################################
######################################################################################################################################################

#########################################################################
### Step 3 - Hierarchical
########################################################################

cat("hierarchical.0.1.2.fss \n")
deriv.hier<-c(0,1,2)
hierarchical.0.1.2.fss<-hierarchical.stepwise.sequential(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, deriv.hier, 
                                                         deriv.analyze[deriv.hier+1], k.grid.in, min.segments.in, seg.sep.in,
                                                         class.method.in, ker.in, seg.weight.in, thresh.in,
                                                         curve.seg.limit.single.in*2, curve.seg.limit.combined.in*2, 
                                                         do.bss=FALSE, density.in, do.par, max.cores, best.sub.max.in)
hierarchical.0.1.2.fss$segment.ids
save.image(output.name)

cat("hierarchical.0.1.2.fss.cv \n")
model.summary.hier.0.1.2.fss<-create.model.list(hierarchical.0.1.2.fss)
hierarchical.0.1.2.fss.cv<-model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, model.summary.hier.0.1.2.fss$derivs, 
                                    indent.deriv[model.summary.hier.0.1.2.fss$derivs+1], model.summary.hier.0.1.2.fss$segment.totals,
                                    model.summary.hier.0.1.2.fss$segments.used, model.summary.hier.0.1.2.fss$npc.param,
                                    folds, trials, folds.list, density.in, do.par, max.cores, large.set=FALSE,
                                    class.method.in, ker.in, seg.weight.in)

save.image(output.name)

###

cat("hierarchical.0.1.2.bss \n")
deriv.hier<-c(0,1,2)
hierarchical.0.1.2.bss<-hierarchical.stepwise.sequential(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, deriv.hier, 
                                                         deriv.analyze[deriv.hier+1], k.grid.in, min.segments.in, seg.sep.in,
                                                         class.method.in, ker.in, seg.weight.in, thresh.in,
                                                         curve.seg.limit.single.in, curve.seg.limit.combined.in, 
                                                         do.bss=TRUE, density.in,do.par, max.cores, best.sub.max.in)
save.image(output.name)

cat("hierarchical.0.1.2.bss.cv \n")
model.summary.hier.0.1.2.bss<-create.model.list(hierarchical.0.1.2.bss)
hierarchical.0.1.2.bss.cv<-model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, model.summary.hier.0.1.2.bss$derivs, 
                                    indent.deriv[model.summary.hier.0.1.2.bss$derivs+1], model.summary.hier.0.1.2.bss$segment.totals,
                                    model.summary.hier.0.1.2.bss$segments.used, model.summary.hier.0.1.2.bss$npc.param,
                                    folds, trials, folds.list, density.in, do.par, max.cores, large.set=FALSE,
                                    class.method.in, ker.in, seg.weight.in)

save.image(output.name)


mean(hierarchical.0.1.2.bss.cv$test.set$accuracy)


#########################################################################
### Step 3 - Combined
########################################################################

cat("combined.global.0.1.2.bss \n")
combined.global.0.1.2.bss<-stepwise.combined(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, deriv.analyze, indent.deriv,
                                                   k.grid.in, density.in, do.par, max.cores, large.set=FALSE,
                                                   class.method.in, ker.in, seg.weight.in, thresh.in,
                                                   min.segments.in, seg.sep.in, do.bss=TRUE, best.sub.max.in, 
                                                   global.npc=TRUE, curve.seg.limit.combined.in, curve.seg.limit.single.in)

save.image(output.name)

cat("combined.global.0.1.2.bss.cv \n")
model.summary.combined.global.0.1.2.bss<-create.model.list(combined.global.0.1.2.bss)
combined.global.0.1.2.bss.cv<-model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, model.summary.combined.global.0.1.2.bss$derivs, 
                                    deriv.indents[model.summary.combined.global.0.1.2.bss$derivs+1], model.summary.combined.global.0.1.2.bss$segment.totals,
                                    model.summary.combined.global.0.1.2.bss$segments.used, model.summary.combined.global.0.1.2.bss$npc.param,
                                    folds, trials, folds.list, density.in, do.par, max.cores, large.set=FALSE,
                                    class.method.in, ker.in, seg.weight.in)

save.image(output.name)

cat("combined.local.0.1.2.bss \n")
combined.local.0.1.2.bss<-stepwise.combined(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, deriv.analyze, indent.deriv,
                                                          k.grid.in, density.in, do.par, max.cores, large.set=FALSE,
                                                          class.method.in, ker.in, seg.weight.in, thresh.in,
                                                          min.segments.in, seg.sep.in, do.bss=TRUE, best.sub.max.in, 
                                                          global.npc=FALSE, curve.seg.limit.combined.in, curve.seg.limit.single.in)

save.image(output.name)

cat("combined.local.0.1.2.bss.cv \n")
model.summary.combined.local.0.1.2.bss<-create.model.list(combined.local.0.1.2.bss)
combined.local.0.1.2.bss.cv<-model.cv(fdata.temp, nbasis.in, fpen.lambda.in, classes.temp, model.summary.combined.local.0.1.2.bss$derivs, 
                                       deriv.indents[model.summary.combined.local.0.1.2.bss$derivs+1], model.summary.combined.local.0.1.2.bss$segment.totals,
                                       model.summary.combined.local.0.1.2.bss$segments.used, model.summary.combined.local.0.1.2.bss$npc.param,
                                       folds, trials, folds.list, density.in, do.par, max.cores, large.set=FALSE,
                                       class.method.in, ker.in, seg.weight.in)

save.image(output.name)



