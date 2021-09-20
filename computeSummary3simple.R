computeSummary3simple<- function(obs.p, perm.p, adj, true.upreg, cell_num, main.describe, 
                           seed.number, cutoff=0.05, perm.size=1000, 
                           auc.values, logFC.df=NULL, statistic="os.t"){
  #setwd("/Users/MelodyJin/Desktop/MAST90108/MDSresearchproject/simulation_plots")
  summ<- as.data.frame(matrix(NA, nrow =ncol(obs.p)*2 , ncol =17))
  colnames(summ)<- c("group","cell.type","num.cells", "test.name", 
                     "num.true.upreg","num.upreg","TP", "precision", 
                     "sensitivity","F1", "overlap","TP.overlap","cutoff",
                     "perm.size","seed.number","AUC","statistic")
  summ[, "group"]<- main.describe
  summ[, "cutoff"]<- cutoff
  summ[, "statistic"]<- statistic
  summ[, "perm.size"]<- perm.size
  summ[, "seed.number"]<- seed.number
  summ[, "cell.type"]<- rep(colnames(obs.p), 2)
  summ[, "test.name"]<- "permutation"
  summ[1:ncol(obs.p), "test.name"]<- rep("limma", ncol(obs.p))
  #summ[ncol(obs.p)+1: ncol(obs.p)*2, "test.name"]<- rep("permutation", ncol(obs.p))
  
  if (statistic == "lfc.p"){
    summ[summ$test.name == "limma", "statistic"]<- "os.t"
  }
  
  for (i in 1: ncol(obs.p)){
    if (substr(statistic,1,2) == "ts"){
      if (is.null(logFC.df)){
        stop("No available logFC data")
      }
      # positive genes 
      limma_upreg_genes =  intersect(row.names(obs.p[which(obs.p[,i] < cutoff),]),
                                     row.names(logFC.df[which(logFC.df[,i] > 0),]))
      perm_upreg_genes =  intersect(row.names(perm.p[which(perm.p[,i] < cutoff),]), 
                                    row.names(logFC.df[which(logFC.df[,i] > 0),]))
      
      # negative genes 
      limma_nonsig_genes =  intersect(row.names(obs.p[which(obs.p[,i] >=cutoff),]),
                                      row.names(logFC.df[which(logFC.df[,i] < 0),]))
      perm_nonsig_genes =  intersect(row.names(perm.p[which(perm.p[,i] >= cutoff),]), 
                                     row.names(logFC.df[which(logFC.df[,i] < 0),]))
      
    }else{
      # positive genes 
      limma_upreg_genes =  row.names(obs.p[which(obs.p[,i] < cutoff),])
      perm_upreg_genes =  row.names(perm.p[which(perm.p[,i] < cutoff),])
      
      # negative genes 
      limma_nonsig_genes =  row.names(obs.p[which(obs.p[,i]>=(1-cutoff)),])
      perm_nonsig_genes =  row.names(perm.p[which(perm.p[,i]>=(1-cutoff)),])
    }
    
    # metrics 
    tp.l<- length(intersect(limma_upreg_genes, true.upreg[[i]]))
    fp.l<- length(limma_upreg_genes)-length(intersect(true.upreg[[i]], limma_upreg_genes))
    tn.l<- length(intersect(setdiff(row.names(obs.p),true.upreg[[i]]),limma_nonsig_genes))
    fn.l<- length(setdiff(true.upreg[[i]], limma_upreg_genes))
    
    tp.p<- length(intersect(perm_upreg_genes, true.upreg[[i]]))
    fp.p<- length(perm_upreg_genes)-length(intersect(true.upreg[[i]], perm_upreg_genes))
    tn.p<- length(intersect(setdiff(row.names(obs.p),true.upreg[[i]]),perm_nonsig_genes))
    fn.p<- length(setdiff(true.upreg[[i]], perm_upreg_genes))
    
    precision.p<- round(tp.p/(tp.p+fp.p),5)
    precision.l<- round(tp.l/(tp.l+fp.l),5)
    # recall=TruePositives / (TruePositives + FalseNegatives)
    sensitivity.p<-  round(tp.p/(tp.p+fn.p),5)
    sensitivity.l<-  round(tp.l/(tp.l+fn.l),5)
    summ[i, "num.cells"]<- cell_num[i]
    summ[i+ncol(obs.p), "num.cells"]<- cell_num[i]
    
    summ[i, "num.true.upreg"]<- length(true.upreg[[i]])
    summ[i+ncol(obs.p), "num.true.upreg"]<- length(true.upreg[[i]])
    
    summ[i, "num.upreg"]<- length(limma_upreg_genes)
    summ[i+ncol(obs.p), "num.upreg"]<- length(perm_upreg_genes)
    
    summ[i, "TP"]<- tp.l
    summ[i+ncol(obs.p), "TP"]<- tp.p
    
    summ[i, "precision"]<- round(tp.l/(tp.l+fp.l),5)
    summ[i+ncol(obs.p), "precision"]<- round(tp.p/(tp.p+fp.p),5)
    
    summ[i, "sensitivity"]<- round(tp.l/(tp.l+fn.l),5)
    summ[i+ncol(obs.p), "sensitivity"]<- round(tp.p/(tp.p+fn.p),5)
    
    summ[i, "F1"]<- round((2 * precision.l * sensitivity.l) / (precision.l + sensitivity.l), 5)
    summ[i+ncol(obs.p), "F1"]<- round((2 * precision.p * sensitivity.p) / (precision.p + sensitivity.p), 5)
    
    summ[i, "overlap"]<- length(intersect(limma_upreg_genes, perm_upreg_genes))
    summ[i+ncol(obs.p), "overlap"]<- length(intersect(limma_upreg_genes, perm_upreg_genes))
    
    summ[i, "TP.overlap"]<- length(intersect(intersect(limma_upreg_genes, true.upreg[[i]]), 
                                             intersect(perm_upreg_genes, true.upreg[[i]])))
    summ[i+ncol(obs.p), "TP.overlap"]<- length(intersect(intersect(limma_upreg_genes, true.upreg[[i]]), 
                                                         intersect(perm_upreg_genes, true.upreg[[i]])))
    
    summ[i, "AUC"]<- auc.values[i]
    summ[i+ncol(obs.p), "AUC"]<- auc.values[i+ncol(obs.p)]
    
  }
  
  return (summ)
}