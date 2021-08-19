plotComparison<-function(sim, perm.size, seed.n, result, plotMA= TRUE, 
                         plotROC = TRUE, plotCumsum = TRUE, groups=2){
  
  lst<- list(permutation = result$t.perm.pval, limma= result$obs.t.pval)
  lst.names<- list(paste("Group1 - ", as.character(table(sim$Group)[1]), " cells",sep=""), 
                   paste("Group2 - ", as.character(table(sim$Group)[2]), " cells",sep=""))
  resp.g1 <- as.vector(ifelse(test=(sim@rowRanges@elementMetadata$DEFacGroup1 > 1), 
                              yes=1, no=0))
  
  resp.g2 <- as.vector(ifelse(test=(sim@rowRanges@elementMetadata$DEFacGroup2 > 1), 
                              yes=1, no=0))
  resp<- list(resp.g1, resp.g2)
  auc.values<- rep(NA, ncol(result$obs.t.pval)*2)
  pdf(paste("seed",seed.number.cp[m],"pdf", sep="."))
  if (plotROC){
  for( i in 1:length(lst)){
    
    pred.perm<-lst[[1]][,i]
    pred.limma<-lst[[2]][,i]
    
    roc_perm <-roc( response=resp[[i]], predictor= pred.perm,plot=FALSE, 
                    legacy.axes=TRUE, percent=FALSE)
    
    
    roc_limma <-roc( response =resp[[i]], predictor=pred.limma, plot=FALSE, 
                     legacy.axes=TRUE, percent=FALSE) 
    
    rocp_df <- as.data.frame(cbind("sensitivity" = roc_perm$sensitivities, 
                                   "FPP" = 1-roc_perm$specificities, 
                                   "AUC" = roc_perm$auc, 
                                   "test.name" =rep("permutation", length(roc_perm$sensitivities))))
    rocl_df <- as.data.frame(cbind("sensitivity" = roc_limma$sensitivities, 
                                   "FPP" = 1- roc_limma$specificities, 
                                   "AUC" = roc_limma$auc, 
                                   "test.name" =rep("limma", length(roc_limma$sensitivities))))
    
    roc_df<- rbind(rocp_df, rocl_df)
    roc_df$sensitivity<- as.numeric(roc_df$sensitivity)
    roc_df$FPP<- as.numeric(roc_df$FPP)
    roc_df$AUC<- as.numeric(roc_df$AUC)
    auc.values[i]<- roc_limma$auc
    auc.values[i+ncol(result$obs.t.pval.adj)]<- roc_perm$auc
    
    plot(ggplot(roc_df) +geom_line(aes(y = sensitivity, x = FPP, color = test.name)) +  
           facet_zoom(x = roc_df$FPP < 0.3)+
           scale_color_manual(labels = c(paste("limma (AUC=", round(unique(roc_df[roc_df$test.name == "limma", "AUC"]),4),")",sep=""), 
                                         paste("permutation (AUC=", round(unique(roc_df[roc_df$test.name == "permutation", "AUC"]),4),")",sep="")),
                              values = c("#F8766D", "#00BFC4"))+
           labs(title=lst.names[[i]], 
                x = "False Positive Rate (1-Specificity)", 
                y = "True Positive Rate (Sensitivity)"))
  }
  #dev.off()
  #pdf(paste("cumsum_plot_seed",seed.number.cp[m],"pdf", sep="."))
  
  }
  
  if (plotCumsum){
  # cumsum plot of adjusuted p-values 
  for (i in 1:ncol(obs.t.pval)){
    pv.obs<-as.data.frame(matrix(0, ncol = 2, nrow = nrow(result$obs.t.pval)))
    pv.perm<- as.data.frame(matrix(0, ncol = 2, nrow = nrow(result$obs.t.pval)))
    pv.obs[, 1]<- result$obs.t.pval[, i]
    pv.perm[, 1]<- result$t.perm.pval[, i]
    pv.obs[, 2]<- 0
    pv.perm[, 2]<- 0
    pv.obs[match(result$true.upreg[[i]],row.names(result$obs.t.pval)),2]<- 1
    pv.perm[match(result$true.upreg[[i]],row.names(result$t.perm.pval)),2]<- 1
    pv.obs<- pv.obs[order(pv.obs[,1]), ]
    pv.perm<- pv.perm[order(pv.perm[,1]), ]
    plot(cumsum(pv.obs[1:500,2]), type="s",main = colnames(result$obs.t.pval)[i], ylab = "cumsum p values")
    lines(cumsum(pv.perm[1:500,2]), type="s",col = "blue")
    abline(a=0, b=1, col="red")
    legend(x = "bottomright",          
           legend = c("limma", "permutation"),
           col = c("black", "blue"),           
           lwd = 2)    
  }
    for (i in 1:ncol(obs.t.pval)){
      pv.obs<-as.data.frame(matrix(0, ncol = 2, nrow = nrow(result$obs.t.pval)))
      pv.perm<- as.data.frame(matrix(0, ncol = 2, nrow = nrow(result$obs.t.pval)))
      pv.obs[, 1]<- result$obs.t.pval.adj[, i]
      pv.perm[, 1]<- result$t.perm.pval.adj[, i]
      pv.obs[, 2]<- 0
      pv.perm[, 2]<- 0
      pv.obs[match(result$true.upreg[[i]],row.names(result$obs.t.pval)),2]<- 1
      pv.perm[match(result$true.upreg[[i]],row.names(result$t.perm.pval)),2]<- 1
      pv.obs<- pv.obs[order(pv.obs[,1]), ]
      pv.perm<- pv.perm[order(pv.perm[,1]), ]
      plot(cumsum(pv.obs[1:500,2]), type="s",main = colnames(obs.t.pval)[i], ylab = "cumsum p values (adjusted)")
      lines(cumsum(pv.perm[1:500,2]), type="s",col = "blue")
      abline(a=0, b=1, col="red")
      legend(x = "bottomright",          
             legend = c("limma", "permutation"),
             col = c("black", "blue"),           
             lwd = 2)    
    }
  }
  # 
  
  
  
  if (plotMA){
  #summary(as.vector(counts(sim[,sim$Group=="Group1"])[match(true.sup.g1, row.names(sim)), ]))
  #summary(as.vector(counts(sim[,sim$Group=="Group2"])[match(true.sup.g1, row.names(sim)), ]))
    
  logcounts.all <- lognormCounts(counts(sim),log=TRUE,prior.count=0.5) 
  all.bct <- factor(sim$Group)
  design <- model.matrix(~0+all.bct)
  colnames(design)[1:(length(levels(all.bct)))] <- levels(all.bct)
  mycont <- matrix(0,ncol=length(levels(all.bct)),nrow=length(levels(all.bct)))
  colnames(mycont)<-levels(all.bct)
  diag(mycont)<-1
  mycont[upper.tri(mycont)]<- -1/(length(levels(all.bct))-1)
  mycont[lower.tri(mycont)]<- -1/(length(levels(all.bct))-1)
  mycont
  # Fill out remaining rows with 0s
  zero.rows <- matrix(0, ncol=length(levels(all.bct)),
                      nrow=(ncol(design)-length(levels(all.bct))))
  test <- rbind(mycont,zero.rows)
  test
  fit.obs <- lmFit(logcounts.all,design)
  fit.contrast <- contrasts.fit(fit.obs,contrasts=test)
  
  fit.cont.eb.obs <- eBayes(fit.contrast,trend=TRUE,robust=TRUE)
  
  dt<-decideTests(fit.cont.eb.obs)
  tp.g1<- row.names(fit.cont.eb.obs) %in% result$true.upreg[[1]]
  tp.g2<- row.names(fit.cont.eb.obs) %in% result$true.upreg[[2]]
  #pdf(paste("plot_MA_seed",seed.number.cp[m],"pdf", sep="."))
  limma::plotMA(fit.cont.eb.obs, coef=1, status =tp.g1, main = "Group1 true DE", 
                legend = "bottomright")
  
  limma::plotMA(fit.cont.eb.obs, coef=2, status =tp.g2, main = "Group2 true DE",
                legend = "bottomright")
  
  limma::plotMA(fit.cont.eb.obs, coef=1, status = dt[,1], main = "Group 1 (limma)")
  
  limma::plotMA(fit.cont.eb.obs, coef=2, status = dt[,2], main = "Group 2 (limma)")
  
  #plotSA(fit.cont.eb.obs)
  
  # plot for permutation
  sig<- row.names(result$t.perm.pval.adj)[result$t.perm.pval.adj$Group1<0.05]
  perm.g1<- row.names(fit.cont.eb.obs) %in% sig #DE genes perm
  limma::plotMA(fit.cont.eb.obs, coef=1, status=perm.g1, main = "Group 1 (permutation)")
  
  sig<- row.names(result$t.perm.pval.adj)[result$t.perm.pval.adj$Group2<0.05]
  perm.g2<- row.names(fit.cont.eb.obs) %in% sig #DE genes perm
  limma::plotMA(fit.cont.eb.obs, coef=2, status=perm.g2, main = "Group 2 (permutation)")
  
  dev.off()
  
  }
  return(auc.values)
  #return(list( raw.p =raw.p, 
  #            adj.p = adj.p ))
  
}