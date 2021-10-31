# cell.type is a factor containing the cell type for each cell 
computeObservation<- function(cm, cell.type, statistic = "os.t", lfc.cutoff=0.1,
                              normalise=FALSE){
  
  #########################
  # Limma DE framework written by Belinda Phipson 
  # normalise the count matrix 
  logcounts.all <- lognormCounts(cm,log=TRUE,prior.count=0.5) 
  all.bct <- factor(cell.type)
  design <- model.matrix(~0+all.bct)
  colnames(design)[1:(length(levels(all.bct)))] <- levels(all.bct)
  mycont <- matrix(0,ncol=length(levels(all.bct)),nrow=length(levels(all.bct)))
  colnames(mycont)<-levels(all.bct)
  diag(mycont)<-1
  mycont[upper.tri(mycont)]<- -1/(length(levels(all.bct))-1)
  mycont[lower.tri(mycont)]<- -1/(length(levels(all.bct))-1)
  # Fill out remaining rows with 0s
  zero.rows <- matrix(0,ncol=length(levels(all.bct)),
                      nrow=(ncol(design)-length(levels(all.bct))))
  test <- rbind(mycont,zero.rows)
  if (normalise==TRUE){
    tic("cyclic.logcounts")
    
    cyclic.logcounts<- limma::normalizeBetweenArrays(logcounts.all, 
                                                     method="cyclicloess")
    toc()
    fit.obs <- limma::lmFit(cyclic.logcounts,design)
    
  }else{
    fit.obs <- limma::lmFit(logcounts.all,design)
  }
  #fit.obs <- limma::lmFit(logcounts.all,design)
  fit.contrast <- limma::contrasts.fit(fit.obs,contrasts=test)
  fit.cont.eb.obs <- limma::eBayes(fit.contrast,trend=TRUE,robust=TRUE)
  #########################
  
  cell.types<-length(levels(all.bct))
  tb<- matrix(nrow = nrow(cm), ncol =cell.types +1)
  colnames(tb)<-c(unlist(strsplit(paste("logFC.g", 1:cell.types,
                                        collapse=" ", sep=""), " ")), 
                  "AveExpr")
  for (j in 1:cell.types){
    tb.ct<- limma::topTable(fit.cont.eb.obs,sort="none",
                            number = nrow(fit.cont.eb.obs), coef=j)
    tb[, j] <-tb.ct[,"logFC"]
  }
  
  tb[, j+1]<-tb.ct[,"AveExpr"]
  row.names(tb)<- row.names(tb.ct)
  tb<-as.data.frame(tb)
  if (statistic == "os.t"){
    # one-sided t statistic 
    obs.stat <- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided t statistic 
    obs.pval<-pt(obs.stat, df=fit.contrast$df.residual, lower.tail=FALSE)
    
  }else if (statistic == "ts.t"){
    obs.stat <- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided t statistic 
    obs.pval<-2*pt(obs.stat, df=fit.contrast$df.residual, lower.tail=FALSE)
    
  }else if (statistic == "os.modt"){
    obs.stat<- fit.cont.eb.obs$t
    obs.pval<- pt(obs.stat, df=fit.cont.eb.obs$df.total, lower.tail=FALSE)
    
  }else if (statistic == "ts.modt"){
    # two-sided moderated t statistic 
    obs.stat<- fit.cont.eb.obs$t
    obs.pval<- fit.cont.eb.obs$p.value
    
  }else if (statistic == "os.treatt"){
    # one-sided treat t statistic 
    obs.stat<- (fit.contrast$coef - lfc.cutoff)/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided treat t statistic 
    obs.pval<-pt(obs.stat, df=fit.contrast$df.residual, lower.tail=FALSE)
    
  }else if (statistic == "ts.treatt"){
    # two-sided treat t statistic 
    treat.result<- treat(fit.cont.eb.obs, lfc=lfc.cutoff, trend=TRUE, robust = TRUE)
    obs.stat<- treat.result$t
    obs.pval<-treat.result$p.value
  }else if (statistic == "lfc.p"){
    
    obs.stat <- fit.contrast$coef*(1-pt(fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma, 
                                     df = fit.contrast$df.residual,lower.tail=FALSE))
    message(paste("Can't compute exact p value for ",statistic, " for limma, use one-sided p value from a t-statistic instead"))
    # one-sided t statistic 
    t.stat <- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided t statistic 
    obs.pval<-pt(t.stat, df=fit.contrast$df.residual, lower.tail=FALSE)
  }else{
    
    message(paste("specified statistic", statistic, " not supported"))
  }
  return (list(obs.stat = obs.stat, 
               obs.pval = obs.pval,  
               fit.cont.eb.obs = fit.cont.eb.obs,
               summary.tb = tb))
  
}