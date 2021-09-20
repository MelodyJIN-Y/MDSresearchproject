computeSummary3<- function(obs.p, perm.p, adj, true.upreg, cell_num, main.describe, 
                           seed.number, cutoff=0.05, perm.size=1000, auc.values,
                           logFC.df=NULL, statistic="os.t"){
  #setwd("/Users/MelodyJin/Desktop/MAST90108/MDSresearchproject/simulation_plots")
  summ<- as.data.frame(matrix(NA, nrow =ncol(obs.p)*2 , ncol =17))
  colnames(summ)<- c("group","cell.type","num.cells", "test.name", "num.true.upreg","num.upreg",
                         "TP", "precision", "sensitivity","F1","overlap","TP.overlap","cutoff","perm.size","seed.number","AUC","statistic")
  summ[, "group"]<- main.describe
  summ[, "cutoff"]<- cutoff
  summ[, "statistic"]<- statistic
  summ[, "perm.size"]<- perm.size
  summ[, "seed.number"]<- seed.number
  summ[, "cell.type"]<- rep(colnames(obs.p), 2)
  summ[, "test.name"]<- "permutation"
  # summ[, "num.cells"]<- rep(as.numeric(summ["#cells", ]), 2)
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
                                      row.names(logFC.df[which(logFC.df[,i] <= 0),]))
      perm_nonsig_genes =  intersect(row.names(perm.p[which(perm.p[,i] >= cutoff),]), 
                                     row.names(logFC.df[which(logFC.df[,i] <= 0),]))
      
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
    
    #summ[i, "FPR"]<- round(fp.l/(fp.l+tn.l), 5)
    #summ[i+ncol(obs.p), "FPR"]<- round(fp.p/(fp.p+tn.p), 5)
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
  if (adj  == ""){
    pdf(paste(main.describe, "sim.pdf",  sep="."))
    
  }else{
    pdf(paste(main.describe, adj,"sim.pdf",  sep="."))
  }
  
  
  for (i in 1:ncol(obs.p)){
    hist(obs.p[,i],
         main = paste(adj, "p-values by limma", "-", colnames(obs.p)[i],"cells"), 
         xlab  =(paste(adj, "p-values")))
    
    hist(perm.p[,i],
         main = paste(adj, "p-values by permutation", "-", colnames(obs.p)[i],"cells"), 
         xlab  =(paste(adj, "p-values")))
    
  }

  df.p<- as.data.frame(rbind(cbind("genes" = row.names(obs.p), "test.names"=rep("limma", nrow(obs.p)), 
                                                "group1" = obs.p[,1], "group2" = obs.p[,2]),
                             cbind("genes" = row.names(perm.p), "test.names" =rep("permutation", nrow(perm.p)),
                                   "group1" = perm.p[,1], "group2" = perm.p[,2])))
  colnames(df.p)<- c("genes", "test.name", "group1", "group2")#,"treu.DE")
  df.p[,3]<- as.numeric(df.p[,3])
  df.p[,4]<- as.numeric(df.p[,4])
  df.pv<- melt(df.p, id = c("test.name","genes"), value.name = "pvalue", variable.name = "group")
  df.pv$true.de<- 0
  df.pv$TP<- 0
  for (i in 1: ncol(obs.p)){
    df.pv[which(df.pv$genes %in% intersect(row.names(obs.p[which(obs.p[,i] < cutoff),]), true.upreg[[i]]) & df.pv$test.name == "limma"& df.pv$group == colnames(df.p)[i+2]), "TP"]<- 1
    df.pv[which(df.pv$genes %in% intersect(row.names(perm.p[which(perm.p[,i] < cutoff),]), true.upreg[[i]]) & df.pv$test.name == "permutation"& df.pv$group == colnames(df.p)[i+2]), "TP"]<- 1
    
  }
  df.pv[which(df.pv$genes %in% true.upreg[[1]] & df.pv$group == "group1"), "true.de"]<- 1
  df.pv[which(df.pv$genes %in% true.upreg[[2]] & df.pv$group == "group2"), "true.de"]<- 1
  # table(true.upreg[[2]] %in% df.pv[df.pv$true.de == 1 & df.pv$group == "group2", "genes"])
  
  plot(ggplot(data=df.pv[df.pv$true.de==1 & df.pv$TP==0, ],
              aes(y=pvalue, col=test.name,x = test.name))+
         geom_boxplot()+
         geom_jitter(position=position_jitter(0.3))+
         facet_grid(~group)+geom_hline(yintercept =0.05, col="black")+
         ylab(paste(adj,"pvalue (true DE)"))+
         theme(legend.position= "none")+
         xlab("")+
         geom_jitter(data = df.pv[which(df.pv$TP==1), 1:6], 
               aes(y =pvalue, x = test.name), 
               position=position_jitter(0.3), color = "darkgreen", fill = "darkgreen"))
  
  plot(ggplot(data=df.pv[df.pv$true.de==1 & df.pv$TP==1, ],
              aes(y=pvalue, col=test.name,
                  x = test.name))+
         geom_boxplot()+
         geom_jitter(position=position_jitter(0.4))+
         facet_grid(~group)+geom_hline(yintercept =1/perm.size, col="black")+
         ylab(paste(adj,"pvalue (TP)"))+
         theme(legend.position= "none")+
         xlab(""))
  

  
  plot(ggplot(df.pv[df.pv$true.de==1, ],
              aes(y=pvalue, col=test.name,x = test.name))+
         geom_violin(draw_quantiles = c(0.5), scale="count")+facet_grid(~group)+
         geom_hline(yintercept =0.05, col="black"))
  
  
  for (i in 1:ncol(obs.p)){
    
    limma_upreg_genes =  row.names(obs.p[which(obs.p[,i]<0.05),])
    perm_upreg_genes =  row.names(perm.p[which(perm.p[,i]<0.05),])
    
    df<- as.data.frame(cbind(obs = obs.p[,i], perm = perm.p[,i]))
    row.names(df)<- row.names(obs.p)
    df$true.DE<- 0
    df[true.upreg[[i]], "true.DE"]<- 1
    
    plot(ggplot(data=df[c(limma_upreg_genes,perm_upreg_genes), ],
                aes(col=factor(true.DE)))+
           geom_jitter(aes(x=df[c(limma_upreg_genes,perm_upreg_genes),"obs"],
                           y=df[c(limma_upreg_genes,perm_upreg_genes), "perm"]),
                       position=position_jitter(0.4))
         +geom_hline(yintercept =0, col="black")+
           ylab(paste(adj,"permutation pvalue"))+
           #theme(legend.position= "none")
      ylab(paste(adj,"permutation pvalue"))+
      xlab(paste(adj,"limma pvalue")))
    #boxplot(perm.p[match( true.upreg[[i]], row.names(perm.p)), ],
    #        data = perm.p, 
    #        main=(paste("up-regulated genes for",colnames(perm.p)[i],"cells by permutation")),
    #        names = colnames(obs.p))
    #abline(h=0.05, lty=2, col = "red", lwd = 2)
    #boxplot(obs.p[match( true.upreg[[i]], row.names(obs.p)), ],
    ##        data = obs.p, 
    #        main=(paste("up-regulated genes for",colnames(perm.p)[i],"cells by limma")),
    #        names = colnames(obs.p))
    #abline(h=0.05, lty=2, col = "red", lwd = 2)
    plot(ggplot(data=df[c(limma_upreg_genes,perm_upreg_genes), ],
                aes(col=factor(true.DE)))+
           geom_jitter(aes(x=df[c(limma_upreg_genes,perm_upreg_genes),"obs"],
                           y=df[c(limma_upreg_genes,perm_upreg_genes), "perm"]),
                       position=position_jitter(0.4))
         +geom_hline(yintercept =0, col="black")+
           ylab(paste(adj,"permutation pvalue"))+
           #theme(legend.position= "none")+
           xlim(0,0.05)+ylim(0,0.05)+
           xlab(paste(adj,"limma pvalue")))

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
     
      
      stripchart(lst.perm,ylab =(paste(adj,"pvalue by permutation")), 
                 data=perm.p,method = "jitter",jitter=0.2,
                 main=(paste("up-regulated genes for",colnames(obs.p)[i],"cells by permutation")),
                 group.names = colnames(obs.p),
                 vertical=TRUE)
      abline(h=0.05, lty=2, col = "red", lwd = 2)
    }

  }
  dev.off()
  return (summ)
}