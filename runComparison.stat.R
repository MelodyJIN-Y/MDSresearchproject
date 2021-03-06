# A slight modified version of runComparison function,
# take sim_perm_arrays as input instead of running real permutation 
runComparison.stat<-function(sim, perm.size, statistic = "os.t", sim_perm_arrays,
                             lfc.cutoff=0.1, n.cores=1, true.upreg, 
                             cyclic.normalise=FALSE){
  message(paste("Statistic = ", statistic))
  obs_res<- computeObs(counts(sim), factor(sim$Group),
                       statistic= statistic, lfc.cutoff=lfc.cutoff)
  obs.stat<- obs_res$obs.stat
  obs.pval<- obs_res$obs.pval
  obs.pval<- as.data.frame(obs.pval)
  cell.types <-levels(factor(sim$Group))
  # permutation stats
  t.perm.array<- sim_perm_arrays
  # two-sided p
  if (substr(statistic,1,2) == "ts"){
    perm.pvals<-apply(expand.grid(x = 1:dim(counts(sim))[1], 
                                  y = 1:dim(t.perm.array)[2]), 1, 
                      function(r) (sum(abs(t.perm.array[r[1],r[2], ]) > abs(obs.stat[r[1],r[2]]))+1)/(perm.size+1) )
    
  }else{
    perm.pvals<-apply(expand.grid(x = 1:dim(counts(sim))[1], 
                                  y = 1:dim(t.perm.array)[2]), 1, 
                      function(r) (sum(t.perm.array[r[1],r[2], ] > obs.stat[r[1],r[2]])+1)/(perm.size+1) )
  }
  
  perm.pval<- matrix(perm.pvals, 
                     nrow=dim(counts(sim))[1], 
                     ncol = dim(t.perm.array)[2], 
                     byrow=FALSE)
  rownames(perm.pval) <- row.names(counts(sim))
  colnames(perm.pval) <- cell.types
  perm.pval.adj<- apply(perm.pval, 2, p.adjust, method="BH")
   

  perm.pval.adj<- as.data.frame(perm.pval.adj)
  rownames(perm.pval.adj) <- row.names(counts(sim))
  colnames(perm.pval.adj) <- cell.types
  
  obs.pval.adj<- apply(obs.pval, 2, p.adjust, method="BH")
  obs.pval.adj <- as.data.frame(obs.pval.adj)
  rownames(obs.pval.adj) <- row.names(counts(sim))
  colnames(obs.pval.adj) <- cell.types
  rownames(obs.pval) <- row.names(counts(sim))
  colnames(obs.pval) <-cell.types
  
  df.list=vector("list", length(cell.types))
  for (ct in 1: length(cell.types)){
    group.df<- obs_res$summary.tb[, c(ct, length(cell.types)+1)] 
    row.names(group.df)<-row.names(obs_res$summary.tb)
    df.list[[ct]] <-group.df
  }
  setNames(df.list, paste("group", 1:ct,".df",sep=""))
  for (i in 1:length(df.list)){
    df.list[[i]]$true.DE<- 0
    df.list[[i]][match(true.upreg[[i]],row.names(obs_res$summary.tb)), "true.DE"]<-1
    df.list[[i]]$limma.raw.pvalue <- obs.pval[,i]
    df.list[[i]]$limma.adj.pvalue <- obs.pval.adj[, i]
    df.list[[i]]$perm.raw.pvalue <- perm.pval[, i]
    df.list[[i]]$perm.adj.pvalue <- perm.pval.adj[, i]
  }
  
  return(list(result.table= df.list , obs.stat= obs.stat, 
              t.perm.array=t.perm.array,
              perm.pval = perm.pval, 
              perm.pval.adj=perm.pval.adj,
              obs.pval=obs.pval,
              obs.pval.adj=obs.pval.adj,
              true.upreg = true.upreg
  ))
}