# A function to run permutation for a given matrix
# Input
# 
pPermutation<- function(cm, groups, perm.size = 1000, statistic = "os.t",  
                 lfc.cutoff=0.1, n.cores=1, seed_value=97){
  
  logcounts.all <- lognormCounts(cm, log=TRUE, prior.count=0.5) 
  all.bct <- factor(groups)
  t.perm.array<- array(0, dim = c(dim(logcounts.all)[1], 
                                  length(unique(all.bct)),perm.size))
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  tic.clearlog()
  tic(paste(perm.size, "permutation total time"))
  doParallel::registerDoParallel(cl = my.cluster)
  set.seed(seed_value)
  sds<- sample(1:200000, size= perm.size, replace = T)
  # for (i in 1:perm.size) {
  t.perm.array<-foreach (i = 1:perm.size) %dopar% {
    # permutate the all.bct only 
    set.seed(sds[i])
    all.bct.shuffled <- sample(all.bct, replace = FALSE, size =length(all.bct))
    
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
    #############################
    
    if (statistic == "os.t"){
      # one-sided t statistic 
      perm.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
      
    }else if (statistic == "ts.t"){
      perm.stat <- fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma
      
    }else if (statistic == "os.modt"){
      perm.stat<- fit.cont.eb$t
      
    }else if (statistic == "ts.modt"){
      # two-sided moderated t statistic 
      perm.stat<- fit.cont.eb$t
    }else if (statistic == "os.treatt"){
      # one-sided treat t statistic 
      perm.stat<- (fit.cont$coef - lfc.cutoff)/ fit.cont$stdev.unscaled / fit.cont$sigma
      
    }else if (statistic == "ts.treatt"){
      # two-sided treat t statistic 
      treat.result<- limma::treat(fit.cont.eb, lfc=lfc.cutoff, trend=TRUE, robust = TRUE)
      perm.stat<- treat.result$t
      
    }else if (statistic == "lfc.p"){
      # logFc * (1-p)
      perm.stat <- fit.cont$coef*(1-pt(fit.cont$coef/ fit.cont$stdev.unscaled / fit.cont$sigma, 
                                       df = fit.cont$df.residual,lower.tail=FALSE))
    }else{
      print(paste("specified statistic", statistic, " not supported"))
    }
    t.perm.array[,,i]<- perm.stat
  }
  toc(log=TRUE)
  perm.array<- simplify2array(t.perm.array)
  parallel::stopCluster(cl = my.cluster)
  
  return (list(ct.names = levels(groups),
               t.perm = perm.array))
}
