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
library(randomForestSRC)

load("myData.RData")



### Code for running lasso in cluster
set.seed(seed)

m = 5

#save results of interest

resnames<-c(paste(true.vars,"chosen",sep="."), paste(true.vars,"surrogates",sep="."),
            "N.tgt.chosen","N.nontgt.chosen", "N.nontgt.chosen.gt.9", "N.nontgt.chosen.gt.8","N.nontgt.chosen.lt.6",
            "N.tgt.sgt.chosen","seed") 

RF.res<-RF_LSO.res<-RF.LSO_val.res<-matrix(rep(NA,m*length(resnames)),nrow=m)
colnames(RF.res)<- resnames
colnames(RF_LSO.res)<-resnames
colnames(RF.LSO_val.res)<-resnames


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
  
  #### Use random forests with resampling of variable importance ####
  XCENwY <- as.data.frame(XCEN[set.ind==0,]) %>%
    mutate(SIMY = as.factor(y[set.ind==0]))
  
  ## grow the forest - request VIMP
  reg.o <- rfsrc(SIMY ~ ., data = XCENwY)
  
  ## subsample - try bootstrap number as 100
  reg.smp.o <- subsample(reg.o, B = 100)
  
  # get 95% LCL
  #subsample.object<-extract.subsample(reg.smp.o)
  #CIs<-subsample.object$ci
  #LCL95pct<-CIs["2.5%",]
  #selected <- names(LCL95pct[LCL95pct>0])
  
  # use Bonferroni adjusted 95% LCL jackknife approach
  subsample.object<-extract.subsample(reg.smp.o)
  SE.Z<-subsample.object$se.Z
  SE.jk<-subsample.object$se.jk.Z
  vmp<-subsample.object$vmp
  
  LCL95pct.BFN.jk<-vmp-qnorm(1-(0.05/1000)/2)*SE.jk #Bonferroni-corrected LCL (jackknife)
  selected <- names(LCL95pct.BFN.jk[LCL95pct.BFN.jk>0])
  
  #Test nonzero covariates in validation set
  sig.res<-data.frame(var=selected)
  rownames(sig.res)<-selected
  sig.res$beta.test<-rep(NA,length(selected))
  sig.res$sig.test<-rep(NA,length(selected))
  for (i in 1:length(selected)) {
    var<-as.character(selected[i])
    mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
    sig.res$beta.test[i]<-summary(mod)$coef[2,1]
    sig.res$sig.test[i]<-summary(mod)$coef[2,4]
    
  }

  val.res<-sig.res[sig.res$sig.test<0.05 & !is.na(sig.res$sig.test),]
  



  #How many of the target variables were identified?
  RF.res[j,"N.tgt.chosen"]<-sum(true.vars %in% val.res$var)
  
  
  n.validated<-dim(val.res)[1]
  
  #And how many non-target variables?
  RF.res[j,"N.nontgt.chosen"]<-n.validated-RF.res[j,"N.tgt.chosen"]
  
  RF.res[j,1:10]<-(true.vars %in% val.res$var)
  

  
  # store the number of non-target variables identified as significant

  other.vars.lasso<-val.res[! val.res$var %in% true.vars,]
  
  if(RF.res[j,"N.nontgt.chosen"]>0) {
  
    Cors.true.sig.lasso<-matrix(rep(NA,(RF.res[j,"N.nontgt.chosen"])*10),nrow=(RF.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(RF.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    RF.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    RF.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    RF.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    RF.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)

  }else {
    RF.res[j,"N.nontgt.chosen.gt.9"]<-0
    RF.res[j,"N.nontgt.chosen.gt.8"]<-0
    RF.res[j,"N.nontgt.chosen.lt.6"]<-0
    RF.res[j,11:20]<-rep(0,10)
  }
  
  RF.res[j,"N.tgt.sgt.chosen"]<-sum(RF.res[j,1:10] | RF.res[j,11:20])
  
  ##### Second step - use LASSO on results of random forests ####
  ##### LASSO with less restrictive (min) penalty ######
  #LASSO
  
 
  cvfit = cv.glmnet(XCEN[set.ind==0,selected], y[set.ind==0],alpha=1,family="binomial")
  lasso.sig.vars.index<-predict(cvfit, s = "lambda.min", type="nonzero")

  
  nonzero.res<-data.frame(
    var=c("",selected[unlist(lasso.sig.vars.index)]),
    coef=coef(cvfit, s = "lambda.min")[which(coef(cvfit, s = "lambda.min") != 0)]
  )
  
  #remove intercept
  nonzero.res<-nonzero.res[nonzero.res$var!="",]
  
  #Test nonzero covariates in validation set
  nonzero.res$beta.test<-rep(NA,dim(nonzero.res)[1])
  nonzero.res$sig.test<-rep(NA,dim(nonzero.res)[1])
  
  for (i in 1:dim(nonzero.res)[1]) {
    var<-as.character(nonzero.res$var[i])
    mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
    nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
    nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
    
  }
  
  #n.nonzero.covariates<-dim(nonzero.res)[1]-1
  
  nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
  
  
  #How many of the target variables were identified?
  RF_LSO.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
  
  
  n.validated<-dim(nonzero.res)[1]
  
  #And how many non-target variables?
  RF_LSO.res[j,"N.nontgt.chosen"]<-n.validated-RF_LSO.res[j,"N.tgt.chosen"]
  
  RF_LSO.res[j,1:10]<-(true.vars %in% nonzero.res$var)
  
  
  
  #  store the number of non-target variables identified as significant

  other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  
  if(RF_LSO.res[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(RF_LSO.res[j,"N.nontgt.chosen"])*10),nrow=(RF_LSO.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(RF_LSO.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    RF_LSO.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    RF_LSO.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    RF_LSO.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    RF_LSO.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    RF_LSO.res[j,"N.nontgt.chosen.gt.9"]<-0
    RF_LSO.res[j,"N.nontgt.chosen.gt.8"]<-0
    RF_LSO.res[j,"N.nontgt.chosen.lt.6"]<-0
    RF_LSO.res[j,11:20]<-rep(0,10)
  }
  RF_LSO.res[j,"N.tgt.sgt.chosen"]<-sum(RF_LSO.res[j,1:10] | RF_LSO.res[j,11:20])
  
  
  
  ##### Use Lasso in validation set
  
  cvfit = cv.glmnet(XCEN[set.ind==1,selected], y[set.ind==1],alpha=1,family="binomial")
  lasso.sig.vars.index<-predict(cvfit, s = "lambda.min", type="nonzero")
  
  
  nonzero.res<-data.frame(
    var=c("",selected[unlist(lasso.sig.vars.index)]),
    coef=coef(cvfit, s = "lambda.min")[which(coef(cvfit, s = "lambda.min") != 0)]
  )
  
  #remove intercept
  nonzero.res<-nonzero.res[nonzero.res$var!="",]
  
  #How many of the target variables were identified?
  RF.LSO_val.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
  
  
  n.validated<-dim(nonzero.res)[1]
  
  #And how many non-target variables?
  RF.LSO_val.res[j,"N.nontgt.chosen"]<-n.validated-RF.LSO_val.res[j,"N.tgt.chosen"]
  
  RF.LSO_val.res[j,1:10]<-(true.vars %in% nonzero.res$var)
  
  
  
  #  store the number of non-target variables identified as significant
  
  other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  
  if(RF.LSO_val.res[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(RF.LSO_val.res[j,"N.nontgt.chosen"])*10),nrow=(RF.LSO_val.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(RF.LSO_val.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    RF.LSO_val.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    RF.LSO_val.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    RF.LSO_val.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    RF.LSO_val.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    RF.LSO_val.res[j,"N.nontgt.chosen.gt.9"]<-0
    RF.LSO_val.res[j,"N.nontgt.chosen.gt.8"]<-0
    RF.LSO_val.res[j,"N.nontgt.chosen.lt.6"]<-0
    RF.LSO_val.res[j,11:20]<-rep(0,10)
  }
  RF.LSO_val.res[j,"N.tgt.sgt.chosen"]<-sum(RF.LSO_val.res[j,1:10] | RF.LSO_val.res[j,11:20])
  
  print(j)
}
end_time = Sys.time()
end_time - start_time 


RF.res[,27]<-seed
RF_LSO.res[,27]<-seed
RF.LSO_val.res[,27]<-seed

# write results to output folder
setwd(resfolder)

write.table(RF.res,"RF_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(RF_LSO.res,"RF_LSO_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(RF.LSO_val.res,"RF_LSO_val_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
