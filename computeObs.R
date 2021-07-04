# cell.type is a factor containing the cell type for each cell 
computeObs<- function(cm, cell.type){
  
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
  
  fit.cont.eb.obs$genes <- ann.keep.all
  
  obs.modtstat<- fit.cont.eb.obs$t
  
  ordinary.t<- fit.contrast$coef/ fit.contrast$stdev.unscaled / fit.contrast$sigma
  obs.tstat<- ordinary.t
  obs.t.pval<-pt(ordinary.t, df=fit.contrast$df.residual, lower.tail=FALSE)
  
  obs.mod.pval<- pt(obs.modtstat, df=fit.cont.eb.obs$df.total, lower.tail=FALSE)
  return (list(obs.tstat = obs.tstat, obs.t.pval = obs.t.pval, 
               obs.modtstat= obs.modtstat, obs.mod.pval= obs.mod.pval, 
               fit.cont.eb.obs = fit.cont.eb.obs))
  
}