# A modified version of pPerm function, will include different statistics for each permutation 
pPermMultipleStatistics<- function(cm, groups, perm.size = 1000, n.cores=1, 
                                   normalise=FALSE){
  if (n.cores>1){
    message(paste("Running", perm.size, "permutation with", n.cores, "cores in parallel"))
  }else{
    message(paste("Running", perm.size, "permutation in sequential"))
  }
  logcounts.all <- lognormCounts(cm, log=TRUE, prior.count=0.5) 
  
  # permutate all.bct 
  all.bct <- factor(groups)
  
  
  # permutate all.bct 
  if (normalise==TRUE){
    tic("cyclic.logcounts")
    
    cyclic.logcounts<- limma::normalizeBetweenArrays(logcounts.all, 
                                                     method="cyclicloess")
    toc()
    
  }
  
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  tic.clearlog()
  tic(paste(perm.size, "permutation total time"))
  doParallel::registerDoParallel(cl = my.cluster)
  join <- function(lst, ...) {
    lapply(seq_along(lst),
           function(i) c(lst[[i]], 
                         lapply(list(...), function(lst2) lst2[[i]])
           )
    )
  }
  result<-list(list(), list(), list())
  result <- foreach(i=1:perm.size, .combine='join', .multicombine=TRUE,
                    .init=list(list(), list(), list())) %dopar% {
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
                                          nrow=(ncol(design)-length(levels(all.bct.shuffled)))
                      )
                      
                      test <- rbind(mycont, zero.rows)
                      #
                      if (normalise==TRUE){
                        fit <- limma::lmFit(cyclic.logcounts,design)
                      }else{
                        fit <- limma::lmFit(logcounts.all,design)
                      }
                      
                      fit.cont <- limma::contrasts.fit(fit, contrasts=test)
                      fit.cont.eb <- limma::eBayes(fit.cont,trend=TRUE,robust=TRUE)
                      table.g1 <- limma::topTable(fit.cont.eb,sort="none",number = nrow(fit.cont.eb), coef=1)
                      table.g2 <- limma::topTable(fit.cont.eb,sort="none",number = nrow(fit.cont.eb), coef=2)
                      
                      # one-sided t statistic 
                      os.t.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
                      
                      #}else if (statistic == "ts.t"){
                      ts.t.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
                      
                      #}else if (statistic == "os.modt"){
                      os.modt.stat<- fit.cont.eb$t
                      
                      #}else if (statistic == "ts.modt"){
                      # two-sided moderated t statistic 
                      ts.modt.stat<- fit.cont.eb$t
                      
                      #}else if (statistic == "os.treatt"){
                      # one-sided treat t statistic 
                     # os.treatt.stat<- (fit.cont$coef - lfc.cutoff)/ fit.cont$stdev.unscaled / fit.cont$sigma
                      
          
                      lfc.p.stat <- fit.cont$coef*(1-pt(fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma, 
                                                        df = fit.cont$df.residual,lower.tail=FALSE))
                      
                      
                      list(os.t.stat, os.modt.stat, lfc.p.stat)
                      
                    }
  toc(log=TRUE)
  #array.lst<-simplify2array(array.lst)
  for (i in 1:length(result)){
    result[[i]]<-simplify2array(result[[i]])
  }
  
  parallel::stopCluster(cl = my.cluster)
  
  #return (list(ct.names = levels(groups), #colnames(fit.cont),
  #             array.lst = array.lst))
  return (list(os.t.stat =result[[1]], 
               os.modt.stat =result[[2]], 
               lfc.p.stat =result[[3]]))
}
