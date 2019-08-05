### Cluster setup

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  
  seed <- 9999999

  setwd("/LocalDir")
  resfolder<-"/LocalDir/results"
}else{
  #for(i in 1:length(args)){
  #  eval(parse(text=args[[i]]))
  #}
  seed<-args[1]
  setwd("/RemoteDir")
  resfolder<-"/RemoteDir/results"
}

### R code
library(glmnet)
library(plyr)
library(ggplot2)
library(dplyr)
library(rje)

load("myData.RData")



### Code for running lasso in cluster
set.seed(seed)

m = 5

#save results of interest

resnames<-c(paste(true.vars,"chosen",sep="."), paste(true.vars,"surrogates",sep="."),
            "N.tgt.chosen","N.nontgt.chosen", "N.nontgt.chosen.gt.9", "N.nontgt.chosen.gt.8","N.nontgt.chosen.lt.6",
            "N.tgt.sgt.chosen","seed","alpha") 

ELNET.1se.res<-ELNET.min.res<-matrix(rep(NA,m*length(resnames)),nrow=m)
colnames(ELNET.1se.res)<-resnames
colnames(ELNET.min.res)<-resnames


#Function used to determine if one of the non-target variables has a high correlation with a target variable
surrogate.chosen<-function(v) {
  max(v)>=0.8
}
alphalist<-seq(from=0.05, to=0.95, by=0.05)

start_time = Sys.time()
for (j in 1:m) {
  
  y<-rep(NA,2000)
  for (i in 1:2000) {
    y[i]<-rbinom(1,1,expit(mu.bin[i]))  
  }

  sim.coef<-rep(NA,1000)
  sim.pval<-rep(NA,1000)
  test.coef<-rep(NA,1000)
  test.pval<-rep(NA,1000)
  # split into discovery and validation datasets
  set.ind<-rbinom(length(y),1,p=0.33333)
  
  x.train<-XCEN[set.ind==0,]
  y.train<-y[set.ind==0]
  
  x.test<-XCEN[set.ind==1,]
  y.test<-y[set.ind==1]
  

  
  #start_time = Sys.time()
  mse_results <- lapply(alphalist, function(x) {
    fit <-  cv.glmnet(x.train, y.train, type.measure="mse", alpha=x,family="binomial")
    yhat <- predict(fit, s=fit$lambda.1se, newx=x.test)
    mse <- mean((y.test - yhat)^2)
    res <- c(x, mse)
    names(res)<-c("alpha","mse")
    return(res)
  })
  #end_time = Sys.time()
  #end_time - start_time
  
  res <- t(as.data.frame(mse_results, check.names = FALSE))
  res <- as.data.frame(res)
  
  min_mse <- res %>% 
    slice(which.min(mse))
  
  #min_mse_each_sim <- rbind(min_mse_each_sim, min_mse)
  
  # now we know the optimal mse in each simulation
  #alpha = min_mse_each_sim[j,1]
  alpha = min_mse$alpha
  ELNET.1se.res[j,"alpha"]<-alpha
  
  #GLMNET with more restrictive penalty
  
  cvfit = cv.glmnet(XCEN[set.ind==0,], y[set.ind==0],alpha=alpha,family="binomial")
  
  
  ##### LASSO with more restrictive (1SE) penalty ######
  lasso.sig.vars.index<-predict(cvfit, s = "lambda.1se", type="nonzero")
  
  
  nonzero.res<-data.frame(
    var=c("",colnames(XCEN)[unlist(lasso.sig.vars.index)]),
    coef=coef(cvfit, s = "lambda.min")[which(coef(cvfit, s = "lambda.1se") != 0)]
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
  ELNET.1se.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
  
  
  n.validated<-dim(nonzero.res)[1]
  
  #And how many non-target variables?
  ELNET.1se.res[j,"N.nontgt.chosen"]<-n.validated-ELNET.1se.res[j,"N.tgt.chosen"]
  
  ELNET.1se.res[j,1:10]<-(true.vars %in% nonzero.res$var)
  

  
  # store the number of non-target variables identified as significant

  other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  
  if(ELNET.1se.res[j,"N.nontgt.chosen"]>0) {
  
    Cors.true.sig.lasso<-matrix(rep(NA,(ELNET.1se.res[j,"N.nontgt.chosen"])*10),nrow=(ELNET.1se.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(ELNET.1se.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    ELNET.1se.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    ELNET.1se.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    ELNET.1se.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    ELNET.1se.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)

  }else {
    ELNET.1se.res[j,"N.nontgt.chosen.gt.9"]<-0
    ELNET.1se.res[j,"N.nontgt.chosen.gt.8"]<-0
    ELNET.1se.res[j,"N.nontgt.chosen.lt.6"]<-0
    ELNET.1se.res[j,11:20]<-rep(0,10)
  }
  
  ELNET.1se.res[j,"N.tgt.sgt.chosen"]<-sum(ELNET.1se.res[j,1:10] | ELNET.1se.res[j,11:20])
  
  
  ##### ELNET with less restrictive (min) penalty ######
  lasso.sig.vars.index<-predict(cvfit, s = "lambda.min", type="nonzero")
  
  
  nonzero.res<-data.frame(
    var=c("",colnames(XCEN)[unlist(lasso.sig.vars.index)]),
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
  ELNET.min.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
  
  
  n.validated<-dim(nonzero.res)[1]
  
  #And how many non-target variables?
  ELNET.min.res[j,"N.nontgt.chosen"]<-n.validated-ELNET.min.res[j,"N.tgt.chosen"]
  
  ELNET.min.res[j,1:10]<-(true.vars %in% nonzero.res$var)
  
  
  
  #  store the number of non-target variables identified as significant

  other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  
  if(ELNET.min.res[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(ELNET.min.res[j,"N.nontgt.chosen"])*10),nrow=(ELNET.min.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(ELNET.min.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    ELNET.min.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    ELNET.min.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    ELNET.min.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    ELNET.min.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    ELNET.min.res[j,"N.nontgt.chosen.gt.9"]<-0
    ELNET.min.res[j,"N.nontgt.chosen.gt.8"]<-0
    ELNET.min.res[j,"N.nontgt.chosen.lt.6"]<-0
    ELNET.min.res[j,11:20]<-rep(0,10)
  }
  ELNET.min.res[j,"N.tgt.sgt.chosen"]<-sum(ELNET.min.res[j,1:10] | ELNET.min.res[j,11:20])
  
  print(j)
}
end_time = Sys.time()
end_time - start_time 

ELNET.1se.res[,27]<-seed
ELNET.min.res[,27]<-seed

# write results to output folder
setwd(resfolder)


write.table(ELNET.1se.res,"ELNET_1se_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(ELNET.min.res,"ELNET_min_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
