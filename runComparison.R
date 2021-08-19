runComparison<-function(sim, perm.size){
  
  obs_res<- computeObs(counts(sim), factor(sim$Group))
  obs.tstat<- obs_res$obs.tstat
  obs.t.pval<- obs_res$obs.t.pval
  obs.t.pval<- as.data.frame(obs.t.pval)

  
  # run permutation
  sim_perm_arrays <- tPerm(counts(sim), sim$Group, perm.size=perm.size)
  
  # permutation stats
  t.perm.array<- sim_perm_arrays$t.perm
  
  # t-stat
  t.perm.pvals<-apply(expand.grid(x = 1:dim(counts(sim))[1], 
                                  y = 1:dim(t.perm.array)[2]), 1, 
                      function(r) (sum(t.perm.array[r[1],r[2], ] > obs.tstat[r[1],r[2]])+1)/(perm.size+1) )
  
  t.perm.pval<- matrix(t.perm.pvals, 
                       nrow=dim(counts(sim))[1], 
                       ncol = dim(t.perm.array)[2], 
                       byrow=FALSE)
  rownames(t.perm.pval) <- row.names(counts(sim))
  colnames(t.perm.pval) <- colnames(sim_perm_arrays$ct.names)
  t.perm.pval.adj<- cbind(p.adjust(t.perm.pval[,1], method="BH"),
                          p.adjust(t.perm.pval[,2], method="BH")
  ) 
  t.perm.pval.adj<- as.data.frame(t.perm.pval.adj)
  rownames(t.perm.pval.adj) <- row.names(counts(sim))
  colnames(t.perm.pval.adj) <- sim_perm_arrays$ct.names
  
  
  obs.t.pval.adj<- cbind(p.adjust(obs.t.pval[,1], method="BH"),
                         p.adjust(obs.t.pval[,2], method="BH")
  )
  obs.t.pval.adj <- as.data.frame(obs.t.pval.adj)
  rownames(obs.t.pval.adj) <- row.names(counts(sim))
  colnames(obs.t.pval.adj) <- sim_perm_arrays$ct.names
  rownames(obs.t.pval) <- row.names(counts(sim))
  colnames(obs.t.pval) <- sim_perm_arrays$ct.names
  
  true.sup.g1<- row.names(sim)[sim@rowRanges@elementMetadata$DEFacGroup1 > 1]
  true.sup.g2<- row.names(sim)[sim@rowRanges@elementMetadata$DEFacGroup2 > 1]
  true.sup.g1<- intersect(true.sup.g1, row.names(sim))
  true.sup.g2<- intersect(true.sup.g2, row.names(sim))
  true.upreg = list(true.sup.g1, true.sup.g2)
  
  group1.df<- obs_res$summary.tb[,c("AveExpr", "logFC.g1")] 
  group2.df<- obs_res$summary.tb[,c("AveExpr", "logFC.g2")] 
  row.names(group1.df)<-row.names(obs_res$summary.tb)
  row.names(group2.df)<-row.names(obs_res$summary.tb)
  df.list<- list(group1.df, group2.df)
  for (i in 1:length(df.list)){
    df.list[[i]]$true.DE<- 0
    df.list[[i]][match(true.upreg[[i]], row.names(group1.df)), "true.DE"]<-1
    df.list[[i]]$limma.raw.pvalue <- obs.t.pval[,i]
    df.list[[i]]$limma.adj.pvalue <- obs.t.pval.adj[, i]
    df.list[[i]]$perm.raw.pvalue <- t.perm.pval[, i]
    df.list[[i]]$perm.adj.pvalue <- t.perm.pval.adj[, i]
    
  }

  return(list(result.table= df.list , obs.tstat= obs.tstat, t.perm.array=t.perm.array,
              t.perm.pval = t.perm.pval, 
              t.perm.pval.adj=t.perm.pval.adj,
              obs.t.pval=obs.t.pval,
              obs.t.pval.adj=obs.t.pval.adj,
              true.upreg = true.upreg
              ))
}