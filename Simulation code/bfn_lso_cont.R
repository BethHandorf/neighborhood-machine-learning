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

  setwd("//LocalDir")
  resfolder<-"//LocalDir//results"
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

# seed<-10003
# seed<-10005
# seed<-10009
# seed<-10012
# seed<-10015
# seed<-10018
# seed<-10020
# seed<-10023
# seed<-10026
# seed<-10027
# seed<-10028
# seed<-10030
 seed<-10032
# seed<-10035
# seed<-10039
# seed<-10042
# seed<-10044
# seed<-10045

### Code for running lasso in cluster
set.seed(seed)

m = 10

#save results of interest

resnames<-c(paste(true.vars,"chosen",sep="."), paste(true.vars,"surrogates",sep="."),
            "N.tgt.chosen","N.nontgt.chosen", "N.nontgt.chosen.gt.9", "N.nontgt.chosen.gt.8","N.nontgt.chosen.lt.6",
            "N.tgt.sgt.chosen","seed") 

BFN.res<-LSO.1se.res<-LSO.min.res<-matrix(rep(NA,m*length(resnames)),nrow=m)
colnames(BFN.res)<- resnames
colnames(LSO.1se.res)<-resnames
colnames(LSO.min.res)<-resnames


#Function used to determine if one of the non-target variables has a high correlation with a target variable
surrogate.chosen<-function(v) {
  max(v)>=0.8
}

start_time = Sys.time()
for (j in 1:m) {
  
  y<-rep(NA,2000)
  for (i in 1:2000) {
    y[i]<-rnorm(1,mu.cont[i],1)  
  }

  sim.coef<-rep(NA,1000)
  sim.pval<-rep(NA,1000)
  test.coef<-rep(NA,1000)
  test.pval<-rep(NA,1000)
  # split into discovery and validation datasets
  set.ind<-rbinom(length(y),1,p=0.33333)
  
  for (i in 1:1000) {
    mod<-glm(y[set.ind==0] ~ XCEN[set.ind==0,][, i])
    sim.coef[i]<-summary(mod)$coef[2,1]
    sim.pval[i]<-summary(mod)$coef[2,4]
    
    mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1,][, i])
    test.coef[i]<-summary(mod)$coef[2,1]
    test.pval[i]<-summary(mod)$coef[2,4]
  }
  
  sim.res<-matrix(cbind(sim.coef,sim.pval,test.coef,test.pval),ncol=4)
  rownames(sim.res)<-colnames(XCEN)
  #What variables have a bonferroni-adjusted significant p-val?
  
  true.vars.res<-sim.res[true.vars,]
  
  true.tmp<-true.vars.res[,2]<(.05/1000) & true.vars.res[,4]<(.05)
  BFN.res[j,1:10]<-true.tmp
  BFN.res[j,"N.tgt.chosen"]<-sum(true.tmp) #total number of "TRUE"
  

  #Other variables
    
  other.vars<-sim.res[! colnames(XCEN) %in% true.vars,]
  other.vars.sig<-other.vars[other.vars[,2]<(.05/1000)& other.vars[,4]<(.05), ,drop=FALSE]
  nsig<-dim(other.vars.sig)[1]
  
  BFN.res[j,"N.nontgt.chosen"]<-nsig
  
  # Store the number of non-target variables identified as significant
  #non.true.vars.res<-sim.res[non.true.vars,]
  #each.non_tgt.BFN.chosen[j,]<-non.true.vars.res[,2]<(.05/1000)
  

  
  if(nsig>0) {
    #What is the correlation between the signifcant variables and each of the 
    #true variables
    Cors.true.sig<-matrix(rep(NA,nsig*10),nrow=nsig)
    rownames(Cors.true.sig)<-rownames(other.vars.sig)
    colnames(Cors.true.sig)<-true.vars
    for (k in 1:10) {
      for (i in 1:nsig){
        Cors.true.sig[i,k]<-cor.abs[rownames(other.vars.sig)[i],true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig<-apply(Cors.true.sig,1,max)
    max.cor.w.sig.varname<-true.vars[apply(Cors.true.sig,1,which.max)]
    
    other.vars.sig<-cbind(other.vars.sig, max.cor.w.sig, max.cor.w.sig.varname)
    
    BFN.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig>=0.9)
    BFN.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig>=0.8)
    BFN.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig<0.6)

    BFN.res[j,11:20]<-apply(Cors.true.sig,2,surrogate.chosen)
    
  
  }else {
    BFN.res[j,"N.nontgt.chosen.gt.9"]<-0
    BFN.res[j,"N.nontgt.chosen.gt.8"]<-0
    BFN.res[j,"N.nontgt.chosen.lt.6"]<-0
    BFN.res[j,11:20]<-rep(0,10)
  }
  
  BFN.res[j,"N.tgt.sgt.chosen"]<-sum(BFN.res[j,1:10] | BFN.res[j,11:20])
  
  #LASSO
  
  cvfit = cv.glmnet(XCEN[set.ind==0,], y[set.ind==0],alpha=1)

  
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
  
  if (dim(nonzero.res)[1]>0) {
  
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var])
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    #n.nonzero.covariates<-dim(nonzero.res)[1]-1
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    LSO.1se.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    LSO.1se.res[j,"N.nontgt.chosen"]<-n.validated-LSO.1se.res[j,"N.tgt.chosen"]
    
    LSO.1se.res[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    # store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  } else{
    LSO.1se.res[j,"N.tgt.chosen"]<-0
    LSO.1se.res[j,"N.nontgt.chosen"]<-0
    LSO.1se.res[j,1:10]<-rep(0,10)  
  }
  
  if(LSO.1se.res[j,"N.nontgt.chosen"]>0) {
  
    Cors.true.sig.lasso<-matrix(rep(NA,(LSO.1se.res[j,"N.nontgt.chosen"])*10),nrow=(LSO.1se.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(LSO.1se.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    LSO.1se.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    LSO.1se.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    LSO.1se.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    LSO.1se.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)

  }else {
    LSO.1se.res[j,"N.nontgt.chosen.gt.9"]<-0
    LSO.1se.res[j,"N.nontgt.chosen.gt.8"]<-0
    LSO.1se.res[j,"N.nontgt.chosen.lt.6"]<-0
    LSO.1se.res[j,11:20]<-rep(0,10)
  }
  
  LSO.1se.res[j,"N.tgt.sgt.chosen"]<-sum(LSO.1se.res[j,1:10] | LSO.1se.res[j,11:20])
  
  
  ##### LASSO with less restrictive (min) penalty ######
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
  
  if (dim(nonzero.res)[1]>0) {
    
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var])
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    #n.nonzero.covariates<-dim(nonzero.res)[1]-1
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    LSO.min.res[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    LSO.min.res[j,"N.nontgt.chosen"]<-n.validated-LSO.min.res[j,"N.tgt.chosen"]
    
    LSO.min.res[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    #  store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  
  } else{
    LSO.min.res[j,"N.tgt.chosen"]<-0
    LSO.min.res[j,"N.nontgt.chosen"]<-0
    LSO.min.res[j,1:10]<-rep(0,10)  
  }
  
    
  if(LSO.min.res[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(LSO.min.res[j,"N.nontgt.chosen"])*10),nrow=(LSO.min.res[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(LSO.min.res[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    LSO.min.res[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    LSO.min.res[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    LSO.min.res[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    LSO.min.res[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    LSO.min.res[j,"N.nontgt.chosen.gt.9"]<-0
    LSO.min.res[j,"N.nontgt.chosen.gt.8"]<-0
    LSO.min.res[j,"N.nontgt.chosen.lt.6"]<-0
    LSO.min.res[j,11:20]<-rep(0,10)
  }
  LSO.min.res[j,"N.tgt.sgt.chosen"]<-sum(LSO.min.res[j,1:10] | LSO.min.res[j,11:20])
  
  print(j)
}
end_time = Sys.time()
end_time - start_time 

BFN.res[,27]<-seed
LSO.1se.res[,27]<-seed
LSO.min.res[,27]<-seed

# write results to output folder
setwd(resfolder)

write.table(BFN.res,"BFN_res_cont_redo.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(LSO.1se.res,"LSO_1se_res_cont_redo.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(LSO.min.res,"LSO_min_res_cont_redo.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
