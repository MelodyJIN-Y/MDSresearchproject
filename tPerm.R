library(limma)
library(tictoc)
tPerm<- function(cm, groups, perm.size = 1000){
  
  logcounts.all <- lognormCounts(cm, log=TRUE, prior.count=0.5) 
  
  # permutate all.bct 
  all.bct <- factor(groups)
  modt.perm.array<- array(0, dim = c(dim(logcounts.all)[1], 
                                     length(unique(all.bct)), perm.size))
  t.perm.array<- array(0, dim = c(dim(logcounts.all)[1], 
                                  length(unique(all.bct)),perm.size))
  
  tic.clearlog()
  tic("1000 perm total")
  
  for (i in 1:perm.size) {
    print(i)
    # permutate the all.bct only 
    all.bct.shuffled <- sample(all.bct, replace = FALSE, size =length(all.bct))
    # to do marker analysis using Limma way
    design <- model.matrix(~0+all.bct.shuffled)
    colnames(design)[1:(length(levels(all.bct.shuffled)))] <- levels(all.bct.shuffled)
    mycont <- matrix(0,
                     ncol=length(levels(all.bct.shuffled)),
                     nrow=length(levels(all.bct.shuffled)))
    
    colnames(mycont)<-levels(all.bct.shuffled)
    diag(mycont)<-1
    mycont[upper.tri(mycont)]<- -1/(length(levels(all.bct.shuffled))-1)
    mycont[lower.tri(mycont)]<- -1/(length(levels(all.bct.shuffled))-1)
    
    # Fill out remaining rows with 0s
    zero.rows <- matrix(0,ncol=length(levels(all.bct.shuffled)),
                        nrow=(ncol(design)-length(levels(all.bct.shuffled))))
    
    test <- rbind(mycont, zero.rows)
    tic("lmFit")
    fit <- lmFit(logcounts.all,design)
    toc(log=TRUE)
    tic("contrasts.fit")
    fit.cont <- contrasts.fit(fit, contrasts=test)
    toc(log=TRUE)
    t.perm.array[,,i]<-fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
    tic("eBayes")
    fit.cont.eb <- eBayes(fit.cont,trend=TRUE,robust=TRUE)
    toc(log=TRUE)
    fit.cont$genes <- ann.keep.all
    # store the statistic for current permutation 
    modt.perm.array[,,i]<- fit.cont.eb$t
    
  }
  toc(log=TRUE)
  
  return (list(ct.names = colnames(fit.cont),t.perm = t.perm.array, modet.perm = modt.perm.array))
}