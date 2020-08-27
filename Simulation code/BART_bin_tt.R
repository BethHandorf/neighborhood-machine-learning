### Cluster setup

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  
  seed <- 9999999

  setwd("R://Lynch//multicollinearity//new sims")
  resfolder<-"R://Lynch//multicollinearity//new sims//results"
}else{
  #for(i in 1:length(args)){
  #  eval(parse(text=args[[i]]))
  #}
  seed<-args[1]
  setwd("/fccc/users/statlab/ehandorf/Lynch/multicollinearity")
  resfolder<-"/fccc/users/statlab/ehandorf/Lynch/multicollinearity/results"
}

### R code
library(glmnet)
library(plyr)
library(ggplot2)
library(dplyr)
library(rje)
#library(randomForestSRC)
#library(BART)
options(java.parameters = "-Xmx5g")
library(bartMachine)


load("myData.RData")



### Code for running lasso in cluster
set.seed(seed)

m = 1

#save results of interest

resnames<-c(paste(true.vars,"chosen",sep="."), paste(true.vars,"surrogates",sep="."),
            "N.tgt.chosen","N.nontgt.chosen", "N.nontgt.chosen.gt.9", "N.nontgt.chosen.gt.8","N.nontgt.chosen.lt.6",
            "N.tgt.sgt.chosen","seed") 

BART.res.local<-BART.res.globalSE<-BART.res.globalMAX<-matrix(rep(NA,m*length(resnames)),nrow=m)
colnames(BART.res.local)<- resnames
colnames(BART.res.globalSE)<-resnames
colnames(BART.res.globalMAX)<-resnames


#Function used to determine if one of the non-target variables has a high correlation with a target variable
surrogate.chosen<-function(v) {
  max(v)>=0.8
}

start_time = Sys.time()
for (j in 1:m) {
  
  y<-rep(NA,2000)
  for (i in 1:2000) {
    y[i]<-rbinom(1,1,expit(mu.bin[i]))  
  }
  
  
  # split into discovery and validation datasets
  set.ind<-rbinom(length(y),1,p=0.33333)
  
  X.train<-XCEN[set.ind==0,]
  Y.train<-as.factor(y[set.ind==0]) #Make the outcome a factor for classification via BART
 
  
  #run the BART model
  bart_machine <- bartMachine(data.frame(X.train), Y.train, mem_cache_for_speed=FALSE)
  
  vs <- var_selection_by_permute(bart_machine,bottom_margin = 10, num_permute_samples = 100, plot=FALSE)
  
  selected.local<-vs$important_vars_local_names
  selected.globalSE<-vs$important_vars_global_se_names
  selected.globalMAX<-vs$important_vars_global_max_names
  
  
  ##### Get summary stats for local criteria
  
  #Test nonzero covariates in validation set
  sig.res<-data.frame(var=selected.local)
  rownames(sig.res)<-selected.local
  sig.res$beta.test<-rep(NA,length(selected.local))
  sig.res$sig.test<-rep(NA,length(selected.local))
  
  if(length(selected.local)>0) {
    for (i in 1:length(selected.local)) {
      var<-as.character(selected.local[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      sig.res$beta.test[i]<-summary(mod)$coef[2,1]
      sig.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
  }
  val.res<-sig.res[sig.res$sig.test<0.05 & !is.na(sig.res$sig.test),]
  

  #How many of the target variables were identified?
  BART.res.local[j,"N.tgt.chosen"]<-sum(true.vars %in% val.res$var)
  
  
  n.validated<-dim(val.res)[1]
  
  #And how many non-target variables?
  BART.res.local[j,"N.nontgt.chosen"]<-n.validated-BART.res.local[j,"N.tgt.chosen"]
  
  BART.res.local[j,1:10]<-(true.vars %in% val.res$var)
  

  
  # store the number of non-target variables identified as significant

  other.vars.lasso<-val.res[! val.res$var %in% true.vars,]
  
  if(BART.res.local[j,"N.nontgt.chosen"]>0) {
  
    Cors.true.sig.lasso<-matrix(rep(NA,(BART.res.local[j,"N.nontgt.chosen"])*10),nrow=(BART.res.local[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(BART.res.local[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    BART.res.local[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    BART.res.local[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    BART.res.local[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    BART.res.local[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)

  }else {
    BART.res.local[j,"N.nontgt.chosen.gt.9"]<-0
    BART.res.local[j,"N.nontgt.chosen.gt.8"]<-0
    BART.res.local[j,"N.nontgt.chosen.lt.6"]<-0
    BART.res.local[j,11:20]<-rep(0,10)
  }
  
  BART.res.local[j,"N.tgt.sgt.chosen"]<-sum(BART.res.local[j,1:10] | BART.res.local[j,11:20])
  
  ##### Get summary stats for globalSE criteria
  
  #Test nonzero covariates in validation set
  sig.res<-data.frame(var=selected.globalSE)
  rownames(sig.res)<-selected.globalSE
  sig.res$beta.test<-rep(NA,length(selected.globalSE))
  sig.res$sig.test<-rep(NA,length(selected.globalSE))
  if(length(selected.globalSE)>0 ){
    for (i in 1:length(selected.globalSE)) {
      var<-as.character(selected.globalSE[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      sig.res$beta.test[i]<-summary(mod)$coef[2,1]
      sig.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
  }
  val.res<-sig.res[sig.res$sig.test<0.05 & !is.na(sig.res$sig.test),]
  
  
  #How many of the target variables were identified?
  BART.res.globalSE[j,"N.tgt.chosen"]<-sum(true.vars %in% val.res$var)
  
  
  n.validated<-dim(val.res)[1]
  
  #And how many non-target variables?
  BART.res.globalSE[j,"N.nontgt.chosen"]<-n.validated-BART.res.globalSE[j,"N.tgt.chosen"]
  
  BART.res.globalSE[j,1:10]<-(true.vars %in% val.res$var)
  
  
  
  # store the number of non-target variables identified as significant
  
  other.vars.lasso<-val.res[! val.res$var %in% true.vars,]
  
  if(BART.res.globalSE[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(BART.res.globalSE[j,"N.nontgt.chosen"])*10),nrow=(BART.res.globalSE[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(BART.res.globalSE[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    BART.res.globalSE[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    BART.res.globalSE[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    BART.res.globalSE[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    BART.res.globalSE[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    BART.res.globalSE[j,"N.nontgt.chosen.gt.9"]<-0
    BART.res.globalSE[j,"N.nontgt.chosen.gt.8"]<-0
    BART.res.globalSE[j,"N.nontgt.chosen.lt.6"]<-0
    BART.res.globalSE[j,11:20]<-rep(0,10)
  }
  
  BART.res.globalSE[j,"N.tgt.sgt.chosen"]<-sum(BART.res.globalSE[j,1:10] | BART.res.globalSE[j,11:20])
  
  
  ##### Get summary stats for globalMax
  
  #Test nonzero covariates in validation set
  sig.res<-data.frame(var=selected.globalMAX)
  rownames(sig.res)<-selected.globalMAX
  sig.res$beta.test<-rep(NA,length(selected.globalMAX))
  sig.res$sig.test<-rep(NA,length(selected.globalMAX))
  if(length(selected.globalMAX)>0){
    for (i in 1:length(selected.globalMAX)) {
      var<-as.character(selected.globalMAX[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      sig.res$beta.test[i]<-summary(mod)$coef[2,1]
      sig.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
  }
  val.res<-sig.res[sig.res$sig.test<0.05 & !is.na(sig.res$sig.test),]
  
  
  #How many of the target variables were identified?
  BART.res.globalMAX[j,"N.tgt.chosen"]<-sum(true.vars %in% val.res$var)
  
  
  n.validated<-dim(val.res)[1]
  
  #And how many non-target variables?
  BART.res.globalMAX[j,"N.nontgt.chosen"]<-n.validated-BART.res.globalMAX[j,"N.tgt.chosen"]
  
  BART.res.globalMAX[j,1:10]<-(true.vars %in% val.res$var)
  
  
  
  # store the number of non-target variables identified as significant
  
  other.vars.lasso<-val.res[! val.res$var %in% true.vars,]
  
  if(BART.res.globalMAX[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(BART.res.globalMAX[j,"N.nontgt.chosen"])*10),nrow=(BART.res.globalMAX[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(BART.res.globalMAX[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    BART.res.globalMAX[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    BART.res.globalMAX[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    BART.res.globalMAX[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    BART.res.globalMAX[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    BART.res.globalMAX[j,"N.nontgt.chosen.gt.9"]<-0
    BART.res.globalMAX[j,"N.nontgt.chosen.gt.8"]<-0
    BART.res.globalMAX[j,"N.nontgt.chosen.lt.6"]<-0
    BART.res.globalMAX[j,11:20]<-rep(0,10)
  }
  
  BART.res.globalMAX[j,"N.tgt.sgt.chosen"]<-sum(BART.res.globalMAX[j,1:10] | BART.res.globalMAX[j,11:20])
  
  print(j)
}
end_time = Sys.time()
end_time - start_time 


BART.res.local[,27]<-seed
BART.res.globalSE[,27]<-seed
BART.res.globalMAX[,27]<-seed

# write results to output folder
setwd(resfolder)

write.table(BART.res.local,"BART_local_bin_p100.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(BART.res.globalSE,"BART_globalSE_bin_p100.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(BART.res.globalMAX,"BART_globalMax_bin_p100.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
