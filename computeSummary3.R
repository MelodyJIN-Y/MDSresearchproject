computeSummary3<- function(obs.p, perm.p, adj, true.upreg, cell_num, main.describe){
  #setwd("/Users/MelodyJin/Desktop/MAST90108/MDSresearchproject/simulation_plots")
  summ<- as.data.frame(matrix(NA, nrow =ncol(obs.p)*2 , ncol =11))
  colnames(summ)<- c("group","cell.type","num.cells", "test.name", "num.true.upreg","num.upreg",
                         "TP", "precision", "sensitivity","FPR","overlap")
  summ[, "group"]<- main.describe
  summ[, "cell.type"]<- rep(colnames(obs.p), 2)
  summ[, "test.name"]<- "permutation"
  # summ[, "num.cells"]<- rep(as.numeric(summ["#cells", ]), 2)
  summ[1:ncol(obs.p), "test.name"]<- rep("limma", ncol(obs.p))
  #summ[ncol(obs.p)+1: ncol(obs.p)*2, "test.name"]<- rep("permutation", ncol(obs.p))
  
  for (i in 1: ncol(obs.p)){
    # positive genes 
    limma_upreg_genes =  row.names(obs.p[which(obs.p[,i] < 0.05),])
    perm_upreg_genes =  row.names(perm.p[which(perm.p[,i] < 0.05),])
    
    # negative genes 
    limma_nonsig_genes =  row.names(obs.p[which(obs.p[,i]>=0.05),])
    perm_nonsig_genes =  row.names(perm.p[which(perm.p[,i]>=0.05),])
    
    # metrics 
    tp.l<- length(intersect(limma_upreg_genes, true.upreg[[i]]))
    fp.l<- length(limma_upreg_genes)-length(intersect(true.upreg[[i]], limma_upreg_genes))
    tn.l<- length(intersect(setdiff(row.names(obs.p),true.upreg[[i]]),limma_nonsig_genes))
    fn.l<- length(setdiff(true.upreg[[i]], limma_upreg_genes))
    
    tp.p<- length(intersect(perm_upreg_genes, true.upreg[[i]]))
    fp.p<- length(perm_upreg_genes)-length(intersect(true.upreg[[i]], perm_upreg_genes))
    tn.p<- length(intersect(setdiff(row.names(obs.p),true.upreg[[i]]),perm_nonsig_genes))
    fn.p<- length(setdiff(true.upreg[[i]], perm_upreg_genes))
    
    summ[i, "num.cells"]<- cell_num[i]
    summ[i+ncol(obs.p), "num.cells"]<- cell_num[i]
    
    summ[i, "num.true.upreg"]<- length(true.upreg[[i]])
    summ[i+ncol(obs.p), "num.true.upreg"]<- length(true.upreg[[i]])
    
    summ[i, "num.upreg"]<- length(limma_upreg_genes)
    summ[i+ncol(obs.p), "num.upreg"]<- length(perm_upreg_genes)
    
    summ[i, "TP"]<- tp.l
    summ[i+ncol(obs.p), "TP"]<- tp.p
    
    summ[i, "precision"]<- 100*round(tp.l/(tp.l+fp.l),5)
    summ[i+ncol(obs.p), "precision"]<- 100*round(tp.p/(tp.p+fp.p),5)
      
    summ[i, "sensitivity"]<- 100*round(tp.l/(tp.l+fn.l),5)
    summ[i+ncol(obs.p), "sensitivity"]<- 100*round(tp.p/(tp.p+fn.p),5)
    
    summ[i, "FPR"]<- 100*round(fp.l/(fp.l+tn.l), 5)
    summ[i+ncol(obs.p), "FPR"]<- 100*round(fp.p/(fp.p+tn.p), 5)
    
    summ[i, "overlap"]<- length(intersect(limma_upreg_genes, perm_upreg_genes))
    summ[i+ncol(obs.p), "overlap"]<- length(intersect(limma_upreg_genes, perm_upreg_genes))
    
  }
  if (adj  == ""){
    pdf(paste(main.describe, "sim.pdf",  sep="."))
    
  }else{
    pdf(paste(main.describe, adj,"sim.pdf",  sep="."))
  }
  
  
  for (i in 1:ncol(obs.p)){
    hist(obs.p[,i],
         main = paste(adj, "p-values by limma", "-", colnames(obs.p)[i],"cells"), 
         xlab  =(paste(adj, "p-values (t)")))
    
    hist(perm.p[,i],
         main = paste(adj, "p-values by permutation", "-", colnames(obs.p)[i],"cells"), 
         xlab  =(paste(adj, "p-values (t)")))
    
  }
  
  for (i in 1:ncol(obs.p)){
    
    limma_upreg_genes =  row.names(obs.p[which(obs.p[,i]<0.05),])
    perm_upreg_genes =  row.names(perm.p[which(perm.p[,i]<0.05),])
    
    plot(perm.p[, i],obs.p[, i],     
         lwd=0.1, pch=".", 
         xlab =(paste(adj, "pvalue by permutation")), 
         ylab=(paste(adj,"pvalue by limma")), 
         main=(paste(colnames(obs.p)[i],"cells")))
    abline(a = 0,b=1, col="blue", lwd = 1.5)
    
    
    plot(perm.p[, i],obs.p[, i],
         pch="*", ylim = c(0,0.05), xlim = c(0,0.05),
         xlab =(paste(adj,"pvalue by permutation")), 
         ylab=(paste(adj,"pvalue by limma ")), col = "orange", 
         main=(paste("significant genes for",colnames(obs.p)[i],"cells")))
    abline(a = 0,b=1, col="blue", lwd = 1.5)
    
    if (length(perm_upreg_genes)>0){
      # plots for observation
      lst.obs<- list("a" = obs.p[match(limma_upreg_genes, row.names(obs.p)), 1],
                     "b" =obs.p[match(limma_upreg_genes, row.names(obs.p)),2 ]
                     #"c" =obs.p[match(limma_upreg_genes, row.names(obs.p)),3 ],
                     #"d" =obs.p[match(limma_upreg_genes, row.names(obs.p)),4 ],
                     #"e" =obs.p[match(limma_upreg_genes, row.names(obs.p)),5 ],
                     #"f" =obs.p[match(limma_upreg_genes, row.names(obs.p)),6 ]
                     )
      
      stripchart(lst.obs,ylab =(paste(adj,"pvalue by limma")), 
                 data = obs.p,method = "jitter",jitter=0.2,
                 main=(paste("up-regulated genes for",colnames(obs.p)[i],"cells by limma")),
                 group.names =colnames(obs.p), 
                 vertical=TRUE)
      abline(h=0.05, lty=2, col = "red", lwd = 2)
      
      # plots for permutation
      lst.perm<- list("a" = perm.p[match(perm_upreg_genes, row.names(perm.p)), 1],
                      "b"=perm.p[match(perm_upreg_genes, row.names(perm.p)),2 ]
                      #"c"=perm.p[match(perm_upreg_genes, row.names(perm.p)),3 ],
                      #"d"=perm.p[match(perm_upreg_genes, row.names(perm.p)),4 ],
                      #"e"=perm.p[match(perm_upreg_genes, row.names(perm.p)),5 ],
                      #"f"=perm.p[match(perm_upreg_genes, row.names(perm.p)),6 ]
                      )
      #boxplot(perm.p[match(perm_upreg_genes, row.names(perm.p)), ],
      #        data = perm.p, 
      #        main=(paste("up-regulated genes for",colnames(perm.p)[i],"cells by limma")),
      #        names = colnames(obs.p))
      #abline(h=0.05, lty=2, col = "red", lwd = 2)
      
      stripchart(lst.perm,ylab =(paste(adj,"pvalue by permutation")), 
                 data=perm.p,method = "jitter",jitter=0.2,
                 main=(paste("up-regulated genes for",colnames(obs.p)[i],"cells by permutation")),
                 group.names = colnames(obs.p),
                 vertical=TRUE)
      abline(h=0.05, lty=2, col = "red", lwd = 2)
    }
    
    
    #dev.off()
  }
  dev.off()
  return (summ)
}