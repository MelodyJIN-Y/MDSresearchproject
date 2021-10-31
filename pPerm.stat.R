# A modified version of pPerm function, will include different statistics for each permutation 
pPerm.stat<- function(cm, groups, perm.size = 1000,lfc.cutoff=0.1, n.cores=1){
  if (n.cores>1){
    message(paste("Running", perm.size, "permutation with", n.cores, "cores in parallel"))
  }else{
    message(paste("Running", perm.size, "permutation in sequential"))
  }
  logcounts.all <- lognormCounts(cm, log=TRUE, prior.count=0.5) 
  
  # permutate all.bct 
  all.bct <- factor(groups)
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
  result<-list(list(), list(), list(), list(), list())
  result <- foreach(i=1:perm.size, .combine='join', .multicombine=TRUE,
                  .init=list(list(), list(), list(), list(),list())) %dopar% {
    # permutate the all.bct only 
    all.bct.shuffled <- sample(all.bct, replace = FALSE, size =length(all.bct))
    
    # to do marker analysis using Limma way
    #########################
    # Limma DE framework written by Belinda Phipson
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
    fit <- limma::lmFit(logcounts.all,design)
    fit.cont <- limma::contrasts.fit(fit, contrasts=test)
    fit.cont.eb <- limma::eBayes(fit.cont,trend=TRUE,robust=TRUE)
    #########################
    table.g1 <- limma::topTable(fit.cont.eb,sort="none",number = nrow(fit.cont.eb), coef=1)
    table.g2 <- limma::topTable(fit.cont.eb,sort="none",number = nrow(fit.cont.eb), coef=2)
    
    # one-sided t statistic 
    os.t.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
    
    ts.t.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
    
    os.modt.stat<- fit.cont.eb$t
    
    # two-sided moderated t statistic 
    ts.modt.stat<- fit.cont.eb$t
    

    os.treatt.stat<- (fit.cont$coef - lfc.cutoff)/ fit.cont$stdev.unscaled / fit.cont$sigma
    

    treat.result<- limma::treat(fit.cont.eb, lfc=lfc.cutoff, trend=TRUE, robust = TRUE)
    ts.treatt.stat<- treat.result$t
    

    lfc.p.stat <- fit.cont$coef*(1-pt(fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma,
                                      df = fit.cont$df.residual,lower.tail=FALSE))
    
    list(os.t.stat, os.modt.stat,os.treatt.stat, ts.modt.stat, lfc.p.stat)
    
  }
  toc(log=TRUE)
 
  for (i in 1:length(result)){
    result[[i]]<-simplify2array(result[[i]])
  }
 
  parallel::stopCluster(cl = my.cluster)
  return (list(os.t.stat =result[[1]], 
               os.modt.stat =result[[2]], 
               os.treatt.stat =result[[3]], 
               ts.modt.stat =result[[4]], 
               lfc.p.stat =result[[5]]))
}
