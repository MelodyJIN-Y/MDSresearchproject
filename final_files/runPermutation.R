runPermuation<-function(cm, group, perm.size, statistic = "os.t", 
                             lfc.cutoff=0.1, n.cores=1, perm.arrays = NULL, 
                        normalise=FALSE){
  message(paste("Statistic = ", statistic))
  #obs_res<- computeObservation(cm, factor(group),statistic= statistic, lfc.cutoff=lfc.cutoff)
  tm1 <- system.time(
    {
      obs_res<- computeObservation(cm, factor(group),statistic= statistic, 
                                   lfc.cutoff=lfc.cutoff, normalise = normalise)
    })
  obs.stat<- obs_res$obs.stat
  obs.pval<- obs_res$obs.pval
  obs.pval<- as.data.frame(obs.pval)
  cell.types <-levels(factor(group))
  if (is.null(perm.arrays) == TRUE){
    if (n.cores>1 ){
      message(paste("Running", perm.size, "permutation with", n.cores, "cores in parallel"))
      message(paste("Estimated waiting time: ", ceiling((perm.size*as.numeric(tm1[3])/n.cores)/60),"minutes"))
    }else{
      message(paste("Running", perm.size, "permutation in sequential"))
      message(paste("Estimated waiting time: ", ceiling(perm.size*as.numeric(tm1[3])/60), "minutes"))
    }
    
    # permutation stats
    perm_stat <- pPermutation(cm, group, perm.size=perm.size, 
                              statistic= statistic,lfc.cutoff=lfc.cutoff, 
                              n.cores=n.cores)
    # permutation stats
    perm.arrays<- perm_stat$t.perm
  }

   # two-sided p
  if (substr(statistic,1,2) == "ts"){
    perm.pvals<-apply(expand.grid(x = 1:dim(cm)[1], 
                                  y = 1:dim(perm.arrays)[2]), 1, 
                      function(r) (sum(abs(perm.arrays[r[1],r[2], ]) > abs(obs.stat[r[1],r[2]]))+1)/(perm.size+1) )
    
  }else{
    perm.pvals<-apply(expand.grid(x = 1:dim(cm)[1], 
                                  y = 1:dim(perm.arrays)[2]), 1, 
                      function(r) (sum(perm.arrays[r[1],r[2], ] > obs.stat[r[1],r[2]])+1)/(perm.size+1) )
  }
  
  perm.pval<- matrix(perm.pvals, 
                     nrow=dim(cm)[1], 
                     ncol = dim(perm.arrays)[2], 
                     byrow=FALSE)
  rownames(perm.pval) <- row.names(cm)
  colnames(perm.pval) <- cell.types
  # multiple testing adjustment 
  perm.pval.adj<- apply(perm.pval, 2, p.adjust, method="BH")
  perm.pval.adj<- as.data.frame(perm.pval.adj)
  rownames(perm.pval.adj) <- row.names(cm)
  colnames(perm.pval.adj) <- cell.types
  
  # multiple testing adjustment 
  obs.pval.adj<- apply(obs.pval, 2, p.adjust, method="BH")
  obs.pval.adj <- as.data.frame(obs.pval.adj)
  rownames(obs.pval.adj) <- row.names(cm)
  colnames(obs.pval.adj) <- cell.types 
  obs.pval <- as.data.frame(obs.pval)
  rownames(obs.pval) <- row.names(cm)
  colnames(obs.pval) <-cell.types
  
  df.list=vector("list", length(cell.types))
  for (ct in 1: length(cell.types)){
    group.df<- obs_res$summary.tb[, c(ct, length(cell.types)+1)] 
    row.names(group.df)<-row.names(obs_res$summary.tb)
    df.list[[ct]] <-group.df
  }
  setNames(df.list, paste("group", 1:ct,".df",sep=""))
  
  for (i in 1:length(df.list)){
    upreg.limma<- row.names(obs.pval.adj)[obs.pval.adj[,i] <0.05]
    upreg.perm<- row.names(perm.pval.adj)[perm.pval.adj[,i] <0.05]
    df.list[[i]]$limma.upreg<- 0
    df.list[[i]]$permutation.upreg<- 0
    df.list[[i]][match(upreg.limma,row.names(obs_res$summary.tb)), "limma.upreg"]<-1
    df.list[[i]][match(upreg.perm,row.names(obs_res$summary.tb)), "permutation.upreg"]<-1
    df.list[[i]]$limma.raw.pvalue <- obs.pval[,i]
    df.list[[i]]$limma.adj.pvalue <- obs.pval.adj[, i]
    df.list[[i]]$perm.raw.pvalue <- perm.pval[, i]
    df.list[[i]]$perm.adj.pvalue <- perm.pval.adj[, i]
  }
  
  return(list(result.table= df.list , obs.stat= obs.stat, 
              perm.arrays=perm.arrays,
              perm.pval = perm.pval, 
              perm.pval.adj=perm.pval.adj,
              obs.pval=obs.pval,
              obs.pval.adj=obs.pval.adj
  ))
}