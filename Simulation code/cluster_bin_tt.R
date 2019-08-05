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
library(pvclust)
library(SGL)
library(Hmisc)

load("myData.RData")



### Code for running lasso in cluster
set.seed(seed)

m = 1

#save results of interest

resnames<-c(paste(true.vars,"chosen",sep="."), paste(true.vars,"surrogates",sep="."),
            "N.tgt.chosen","N.nontgt.chosen", "N.nontgt.chosen.gt.9", "N.nontgt.chosen.gt.8","N.nontgt.chosen.lt.6",
            "N.tgt.sgt.chosen","seed","nclst") 

gl.corr<-sgl.corr<-gl.boot<-sgl.boot<-matrix(rep(NA,m*length(resnames)),nrow=m)
#colnames(boot.pca)<-resnames
colnames(gl.corr)<-resnames
colnames(sgl.corr)<-resnames
colnames(gl.boot)<-resnames
colnames(sgl.boot)<-resnames


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
  
  #### #hierarchcal clustering with cutoff at rho=0.8
  #Define "distance" as 1- absolute value of correlation
  corr.spearman.train=cor(x.train,method="spearman")
  
  cor.abs.train = abs(corr.spearman.train)
  
  clst.obj<-hclust(as.dist(1-cor.abs.train),method="complete")
  
  #clst.obj<-hclust(as.dist(1-cor.abs),method="ward.D2") - DONT use
  plot(clst.obj,labels = FALSE)
  clusters<-cutree(clst.obj ,h=0.2)  #Use 0.2 - somwhere invisibly this correlation is being subtracted from 1
  
  #number of clusters
  gl.corr[j,"nclst"]<-length(unique(clusters))
  sgl.corr[j,"nclst"]<-length(unique(clusters))
  
  #Sort columns of XCEN by cluster number
  XCEN.clord<-x.train[,order(clusters)]
  clusters.ord<-clusters[order(clusters)]
  
  ###########################################
  ##### clustering using pvclust function
  

  res.pv2 <- pvclust(XCEN[set.ind==0,], method.dist="abscor", 
                     method.hclust="complete", nboot = 100,
                     parallel=FALSE)
  clusters2 <- pvpick(res.pv2)
 # clusters2
  

  #Variables in a significant cluster
  cluster_vars_list <- unlist(clusters2[[1]])
  
  #Variables not in a significant cluster
  unclst_var_list<-colnames(x.train)[!(colnames(x.train) %in% cluster_vars_list)]
  all_unclust<-x.train[,unclst_var_list]
  
  
  #Generate the 1st PC for each cluster
  clusters_all<-clusters2[[1]]
  allpc <- as.data.frame(matrix(nrow = dim(XCEN[set.ind==0,])[1]))
  
  clusters.boot<-numeric(0)
  
  #file with variable names and bootstrap based cluster membership
  boot.clust.memb<-data.frame(var=colnames(x.train),clust=rep(NA,1000))
  
  for (i in 1:length(clusters_all)) {
    c <- unlist(clusters_all[i])
    c.xcen <- XCEN[set.ind==0,][,c]
    
    pc <-prcomp(c.xcen, scale=TRUE)
    pred <- predict(pc)
    
    firstpc <- as.data.frame(pred[,1])
    allpc <- cbind(allpc,firstpc)
    
    #Save the cluster number with the name
    boot.clust.memb$clust[boot.clust.memb$var %in% c]<-i
  }
  
  allpc <- allpc[,2:(length(clusters_all)+1)]
  colnames(allpc) <- 1:length(clusters_all)
  
  #matrix with PCs for clusters, underlying values for unclustered variables
  clusters.boot.pc<-cbind(allpc,all_unclust)
  
  #Assign cluster numbers to singletons
  boot.clust.memb<-boot.clust.memb[order(boot.clust.memb$clust),]
  seq.tmp<-seq(from=length(clusters_all)+1, to=length(clusters_all)+sum(is.na(boot.clust.memb$clust)),by=1)
  boot.clust.memb$clust[is.na(boot.clust.memb$clust)]<-seq.tmp
  
  XCEN.clord.boot<-x.train[,order(boot.clust.memb$clust)]
  
  
  #Save the number of clusters
 # boot.pca[j,"nclst"]<-dim(clusters.boot)[2]
  gl.boot[j,"nclst"]<-dim(clusters.boot.pc)[2]
  sgl.boot[j,"nclst"]<-dim(clusters.boot.pc)[2]
  
  ################################
  #Clusters are made.  Now fit the algorithms
  ################################
  
  
  #####Group Lasso, clustering at corr>0.8
  
  #Unstandardized group lasso penalty
  cv.GL<-cvSGL(data=list(x=XCEN.clord, y=y.train), type="logit",index=clusters.ord,  
               alpha=0, standardize=FALSE)
  cv.GL.betas<-cv.GL$`fit`$`beta`[,which.min(cv.GL$lldiff) ]
  #names(cv.GL.betas)<-colnames(XCEN.clord)
  
  allres<-data.frame(var=colnames(XCEN.clord),beta=cv.GL.betas,cluster=clusters.ord)
  nonzero.res<-allres[allres$beta!=0,]
  
  #Test nonzero covariates in validation set
   nonzero.res$beta.test<-rep(NA,dim(nonzero.res)[1])
   nonzero.res$sig.test<-rep(NA,dim(nonzero.res)[1])
   
  if (dim(nonzero.res)[1]>0) {
  
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    gl.corr[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    gl.corr[j,"N.nontgt.chosen"]<-n.validated-gl.corr[j,"N.tgt.chosen"]
    
    gl.corr[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    # store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  } else{
    gl.corr[j,"N.tgt.chosen"]<-0
    gl.corr[j,"N.nontgt.chosen"]<-0
    gl.corr[j,1:10]<-rep(0,10)  
  }
  
  if(gl.corr[j,"N.nontgt.chosen"]>0) {
  
    Cors.true.sig.lasso<-matrix(rep(NA,(gl.corr[j,"N.nontgt.chosen"])*10),nrow=(gl.corr[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(gl.corr[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    gl.corr[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    gl.corr[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    gl.corr[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    gl.corr[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)

  }else {
    gl.corr[j,"N.nontgt.chosen.gt.9"]<-0
    gl.corr[j,"N.nontgt.chosen.gt.8"]<-0
    gl.corr[j,"N.nontgt.chosen.lt.6"]<-0
    gl.corr[j,11:20]<-rep(0,10)
  }
  
  gl.corr[j,"N.tgt.sgt.chosen"]<-sum(gl.corr[j,1:10] | gl.corr[j,11:20])
  
  
 
  #####Sparse group Lasso, clustering at corr>0.8
  
  #Unstandardized group lasso penalty
  cv.SGL<-cvSGL(data=list(x=XCEN.clord, y=y.train), type="logit",index=clusters.ord,  
               alpha=0.95, standardize=FALSE)
  cv.SGL.betas<-cv.SGL$`fit`$`beta`[,which.min(cv.SGL$lldiff) ]
  #names(cv.GL.betas)<-colnames(XCEN.clord)
  
  allres<-data.frame(var=colnames(XCEN.clord),beta=cv.SGL.betas,cluster=clusters.ord)
  nonzero.res<-allres[allres$beta!=0,]
  
  #Test nonzero covariates in validation set
  nonzero.res$beta.test<-rep(NA,dim(nonzero.res)[1])
  nonzero.res$sig.test<-rep(NA,dim(nonzero.res)[1])
  
  if (dim(nonzero.res)[1]>0) {
    
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    sgl.corr[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    sgl.corr[j,"N.nontgt.chosen"]<-n.validated-sgl.corr[j,"N.tgt.chosen"]
    
    sgl.corr[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    # store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  } else{
    sgl.corr[j,"N.tgt.chosen"]<-0
    sgl.corr[j,"N.nontgt.chosen"]<-0
    sgl.corr[j,1:10]<-rep(0,10)  
  }
  
  if(sgl.corr[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(sgl.corr[j,"N.nontgt.chosen"])*10),nrow=(sgl.corr[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(sgl.corr[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    sgl.corr[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    sgl.corr[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    sgl.corr[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    sgl.corr[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    sgl.corr[j,"N.nontgt.chosen.gt.9"]<-0
    sgl.corr[j,"N.nontgt.chosen.gt.8"]<-0
    sgl.corr[j,"N.nontgt.chosen.lt.6"]<-0
    sgl.corr[j,11:20]<-rep(0,10)
  }
  
  sgl.corr[j,"N.tgt.sgt.chosen"]<-sum(sgl.corr[j,1:10] | sgl.corr[j,11:20])
  
  #####Group Lasso, bootstrap-based significant clusters only
  
  #Unstandardized group lasso penalty
  cv.GL.boot<-cvSGL(data=list(x=XCEN.clord.boot, y=y.train), type="logit",index=boot.clust.memb$clust,  
               alpha=0, standardize=FALSE)
  cv.GL.boot.betas<-cv.GL.boot$`fit`$`beta`[,which.min(cv.GL.boot$lldiff) ]
  #names(cv.GL.betas)<-colnames(XCEN.clord)
  
  allres<-data.frame(var=colnames(XCEN.clord.boot),beta=cv.GL.boot.betas,cluster=boot.clust.memb$clust)
  nonzero.res<-allres[allres$beta!=0,]
  
  #Test nonzero covariates in validation set
  nonzero.res$beta.test<-rep(NA,dim(nonzero.res)[1])
  nonzero.res$sig.test<-rep(NA,dim(nonzero.res)[1])
  
  if (dim(nonzero.res)[1]>0) {
    
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    gl.boot[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    gl.boot[j,"N.nontgt.chosen"]<-n.validated-gl.boot[j,"N.tgt.chosen"]
    
    gl.boot[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    # store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  } else{
    gl.boot[j,"N.tgt.chosen"]<-0
    gl.boot[j,"N.nontgt.chosen"]<-0
    gl.boot[j,1:10]<-rep(0,10)  
  }
  
  if(gl.boot[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(gl.boot[j,"N.nontgt.chosen"])*10),nrow=(gl.boot[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(gl.boot[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    gl.boot[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    gl.boot[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    gl.boot[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    gl.boot[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    gl.boot[j,"N.nontgt.chosen.gt.9"]<-0
    gl.boot[j,"N.nontgt.chosen.gt.8"]<-0
    gl.boot[j,"N.nontgt.chosen.lt.6"]<-0
    gl.boot[j,11:20]<-rep(0,10)
  }
  
  gl.boot[j,"N.tgt.sgt.chosen"]<-sum(gl.boot[j,1:10] | gl.boot[j,11:20])
  
  
  ##### Sparse Group Lasso, bootstrap-based significant clusters only
  
  #Unstandardized group lasso penalty
  cv.SGL.boot<-cvSGL(data=list(x=XCEN.clord.boot, y=y.train), type="logit",index=boot.clust.memb$clust,  
                    alpha=0.95, standardize=FALSE)
  cv.SGL.boot.betas<-cv.SGL.boot$`fit`$`beta`[,which.min(cv.SGL.boot$lldiff) ]
  #names(cv.GL.betas)<-colnames(XCEN.clord)
  
  allres<-data.frame(var=colnames(XCEN.clord.boot),beta=cv.SGL.boot.betas,cluster=boot.clust.memb$clust)
  nonzero.res<-allres[allres$beta!=0,]
  
  #Test nonzero covariates in validation set
  nonzero.res$beta.test<-rep(NA,dim(nonzero.res)[1])
  nonzero.res$sig.test<-rep(NA,dim(nonzero.res)[1])
  
  if (dim(nonzero.res)[1]>0) {
    
    for (i in 1:dim(nonzero.res)[1]) {
      var<-as.character(nonzero.res$var[i])
      mod<-glm(y[set.ind==1] ~ XCEN[set.ind==1, var],family="binomial")
      nonzero.res$beta.test[i]<-summary(mod)$coef[2,1]
      nonzero.res$sig.test[i]<-summary(mod)$coef[2,4]
      
    }
    
    nonzero.res<-nonzero.res[nonzero.res$sig.test<0.05 & !is.na(nonzero.res$sig.test),]
    
    
    #How many of the target variables were identified?
    sgl.boot[j,"N.tgt.chosen"]<-sum(true.vars %in% nonzero.res$var)
    
    
    n.validated<-dim(nonzero.res)[1]
    
    #And how many non-target variables?
    sgl.boot[j,"N.nontgt.chosen"]<-n.validated-sgl.boot[j,"N.tgt.chosen"]
    
    sgl.boot[j,1:10]<-(true.vars %in% nonzero.res$var)
    
    
    
    # store the number of non-target variables identified as significant
    
    other.vars.lasso<-nonzero.res[! nonzero.res$var %in% true.vars,]
  } else{
    sgl.boot[j,"N.tgt.chosen"]<-0
    sgl.boot[j,"N.nontgt.chosen"]<-0
    sgl.boot[j,1:10]<-rep(0,10)  
  }
  
  if(sgl.boot[j,"N.nontgt.chosen"]>0) {
    
    Cors.true.sig.lasso<-matrix(rep(NA,(sgl.boot[j,"N.nontgt.chosen"])*10),nrow=(sgl.boot[j,"N.nontgt.chosen"]))
    rownames(Cors.true.sig.lasso)<- as.character(other.vars.lasso$var)
    colnames(Cors.true.sig.lasso)<-true.vars
    for (k in 1:10) {
      for (i in 1:(sgl.boot[j,"N.nontgt.chosen"])){
        Cors.true.sig.lasso[i,k]<-cor.abs[as.character(other.vars.lasso$var[i]),true.vars[k]]
      }
    }
    
    #What's the maximum correlation, and with which variables?
    max.cor.w.sig.lasso<-apply(Cors.true.sig.lasso,1,max)
    
    sgl.boot[j,"N.nontgt.chosen.gt.9"]<-sum(max.cor.w.sig.lasso>=0.9)
    sgl.boot[j,"N.nontgt.chosen.gt.8"]<-sum(max.cor.w.sig.lasso>=0.8)
    sgl.boot[j,"N.nontgt.chosen.lt.6"]<-sum(max.cor.w.sig.lasso<.6)
    
    
    surrogate.chosen<-function(v) {
      max(v)>=0.8
    }
    
    sgl.boot[j,11:20]<-apply(Cors.true.sig.lasso,2,surrogate.chosen)
    
  }else {
    sgl.boot[j,"N.nontgt.chosen.gt.9"]<-0
    sgl.boot[j,"N.nontgt.chosen.gt.8"]<-0
    sgl.boot[j,"N.nontgt.chosen.lt.6"]<-0
    sgl.boot[j,11:20]<-rep(0,10)
  }
  
  sgl.boot[j,"N.tgt.sgt.chosen"]<-sum(sgl.boot[j,1:10] | sgl.boot[j,11:20])
  
  print(j)
}
end_time = Sys.time()
end_time - start_time 



gl.corr[,27]<-seed
sgl.corr[,27]<-seed
gl.boot[,27]<-seed
sgl.boot[,27]<-seed
# write results to output folder
setwd(resfolder)


write.table(gl.corr,"gl_corr_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(sgl.corr,"sgl_corr_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(gl.boot,"gl_boot_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)
write.table(sgl.boot,"sgl_boot_res_bin.csv", col.names =FALSE,append=TRUE, sep = ",", row.names = FALSE,quote=FALSE)


