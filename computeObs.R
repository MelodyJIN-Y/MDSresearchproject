# cell.type is a factor containing the cell type for each cell 
computeObs<- function(cm, cell.type, statistic = "t", lfc.cutoff=0.1){
  
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
  
  fit.obs <- lmFit(logcounts.all,design)
  fit.contrast <- contrasts.fit(fit.obs,contrasts=test)
  
  fit.cont.eb.obs <- eBayes(fit.contrast,trend=TRUE,robust=TRUE)
  # fit.cont.eb.obs$genes <- ann.keep.all
  
  if (statistic == "os.t"){
    # one-sided t statistic 
    obs.stat <- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided t statistic 
    obs.pval<-pt(ordinary.t, df=fit.contrast$df.residual, lower.tail=FALSE)
  }else if (statistic == "ts.t"){
    obs.stat <- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
    # p value for one-sided t statistic 
    obs.pval<-2*pt(ordinary.t, df=fit.contrast$df.residual, lower.tail=FALSE)
    
    
    
  }else if (statistic == "os.modt"){
    obs.stat<- fit.cont.eb.obs$t
    obs.mod.pval<- pt(obs.modtstat, df=fit.cont.eb.obs$df.total, lower.tail=FALSE)
    
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
    obs.stat<- treat(fit.cont.eb.obs, lfc=lfc.cutoff, trend=TRUE, robust = TRUE)
    obs.pval<-obs.stat$p.value
    
    
    
  }else if (statistic == "lfc.p"){
    # logFc * (1-p)
    new.t <- fit.contrast$coef*(1-pt(fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma, 
                                     df = fit.contrast$df.residual,lower.tail=FALSE))
  }
  

  table.g1 <- topTable(fit.cont.eb.obs,sort="none",number = nrow(fit.cont.eb.obs), coef=1)
  table.g2 <- topTable(fit.cont.eb.obs,sort="none",number = nrow(fit.cont.eb.obs), coef=2)
  
  tb<- as.data.frame(cbind(logFC.g1 = table.g1[, c("logFC")],
                           logFC.g2 = table.g2[, c("logFC")], AveExpr = table.g2[,"AveExpr"]))
  row.names(tb)<- row.names(table.g1)
  
  
  return (list(obs.stat = obs.stat, 
               obs.pval = obs.pval,  
               fit.cont.eb.obs = fit.cont.eb.obs,
               summary.tb = tb))
  
}