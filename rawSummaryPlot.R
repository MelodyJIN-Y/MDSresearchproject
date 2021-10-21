rawSummaryPlots<- function(plt.title, result, markers, obs, statistic){
  pdf(paste(plt.title, "pdf", sep="."))
  auc.values<- as.data.frame(matrix(0, nrow=2, ncol=ncol(result$obs.pval)+1))
  colnames(auc.values)<- c("test", colnames(result$obs.pval))
  auc.values$test<- c("limma","permutation")
  for (ct in 1:ncol(result$obs.pval)){
    curr.ct.name<- colnames(result$obs.pval)[ct]
    if (curr.ct.name %in% names(markers)){
      # cumsum plot
      curr.ct.name<- colnames(result$obs.pval)[ct]
      pv.obs<-as.data.frame(matrix(0, ncol = 4, nrow = nrow(result$obs.pval)))
      pv.perm<- as.data.frame(matrix(0, ncol = 4, nrow = nrow(result$obs.pval)))
      pv.obs[, 1]<- result$obs.pval[, ct]
      pv.perm[, 1]<- result$perm.pval[, ct]
      pv.obs[, 2]<- 0 # limma TP
      pv.perm[, 2]<- 0 # permutation TP
      pv.obs[, 3]<- 0 # limma FP
      pv.perm[, 3]<- 0 # permutation FP 
      
      pv.obs[, 4]<- 0 # limma t statistic
      pv.perm[, 4]<- 0 # permutation t statistic
      
      pv.obs[match(intersect(markers[[curr.ct.name]],row.names(result$obs.pval)), 
                   row.names(result$obs.pval)) ,2]<- 1
      pv.perm[match(intersect(markers[[curr.ct.name]],row.names(result$perm.pval)), 
                    row.names(result$obs.pval)),2]<- 1
      pv.obs[match(setdiff(row.names(result$obs.pval)[result$obs.pval[,ct]< 0.05], 
                           markers[[curr.ct.name]]),
                   row.names(result$perm.pval)),3]<- 1
      pv.perm[match(setdiff(row.names(result$perm.pval)[result$perm.pval[,ct]< 0.05],
                            markers[[curr.ct.name]]),
                    row.names(result$perm.pval)),3]<- 1
      
      pv.perm[match(row.names(obs$obs.stat),
                    row.names(result$perm.pval)),4]<- obs$obs.stat[match(row.names(obs$obs.stat),
                                                                         row.names(result$perm.pval)), curr.ct.name]
      
      pv.obs[match(row.names(obs$obs.stat),
                   row.names(result$obs.pval)),4]<- obs$obs.stat[match(row.names(obs$obs.stat),
                                                                       row.names(result$obs.pval)), curr.ct.name]
      # sort by p-value, and then by TP
      pv.obs<- pv.obs[order(pv.obs[,1], -pv.obs[,4]), ]
      pv.perm<- pv.perm[order(pv.perm[,1], -pv.perm[,4]), ]
      plot(cumsum(pv.obs[1:500,2]), type="s",
           main = paste(curr.ct.name,"cells",sep=" "), 
           ylab = "cumulative number of TP",lwd=2)
      lines(cumsum(pv.perm[1:500,2]), type="s",col = "blue", lwd=2)
      abline(a=0, b=1, col="red")
      legend(x = "bottomright",          
             legend = c("limma", "permutation"),
             col = c("black", "blue"),           
             lwd = 2)    
      #plot(cumsum(pv.obs[1:10000,2]), type="s",
      #     main = paste(curr.ct.name,"cells",sep=" "), 
      #     ylab = "cumulative number of TP",lwd=2)
      #lines(cumsum(pv.perm[1:10000,2]), type="s",col = "blue", lwd=2)
      #abline(a=0, b=1, col="red")
      #legend(x = "bottomright",          
      #       legend = c("limma", "permutation"),
      #       col = c("black", "blue"),           
      #       lwd = 2)    
      
      pv.obs<- pv.obs[order(pv.obs[,1], -pv.obs[,4]), ]
      pv.perm<- pv.perm[order(pv.perm[,1], -pv.perm[,4]), ]
      plot(cumsum(pv.obs[1:200,3]), type="s",
           main = paste(curr.ct.name,"cells",sep=" "), 
           ylab = "cumulative number of FP",lwd=2)
      lines(cumsum(pv.perm[1:200,3]), type="s",col = "blue",lwd=2)
      #abline(a=0, b=1, col="red")
      legend(x = "topleft",          
             legend = c("limma", "permutation"),
             col = c("black", "blue"),           
             lwd = 2)  
      ################
      
      
      #AUC plot
      resp=list()
      
      resp.g <- as.vector(ifelse(test=(row.names(result$obs.pval) %in% markers[[curr.ct.name]]),  
                                 yes=1, no=0))
      resp[[curr.ct.name]]<-resp.g
      
      pred.perm<-result$perm.pval[,curr.ct.name]
      pred.limma<-result$obs.pval[,curr.ct.name]
      
      roc_perm <-roc( response=resp[[curr.ct.name]], predictor= pred.perm,plot=FALSE, 
                      legacy.axes=TRUE, percent=FALSE)
      
      
      roc_limma <-roc( response =resp[[curr.ct.name]], predictor=pred.limma, plot=FALSE, 
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
      auc.values[auc.values$test=="limma",curr.ct.name]<- roc_limma$auc
      auc.values[auc.values$test=="permutation",curr.ct.name]<- roc_perm$auc
      
      plot(ggplot(roc_df) +geom_line(aes(y = sensitivity, x = FPP, color = test.name)) +  
             facet_zoom(x = roc_df$FPP < 0.3)+
             scale_color_manual(labels = c(paste("limma (AUC=", round(unique(roc_df[roc_df$test.name == "limma", "AUC"]),4),")",sep=""), 
                                           paste("permutation (AUC=", round(unique(roc_df[roc_df$test.name == "permutation", "AUC"]),4),")",sep="")),
                                values = c("#F8766D", "#00BFC4"))+
             labs(title=curr.ct.name, 
                  x = "False Positive Rate (1-Specificity)", 
                  y = "True Positive Rate (Sensitivity)"))
      
      
      
      curr_df<- data.frame(cbind("Amean"=obs$fit.cont.eb.obs$Amean, 
                                 "coef"=obs$fit.cont.eb.obs$coefficients[,ct],
                                 "limma.pvalue"=as.numeric(result$obs.pval[,ct]),
                                 "perm.pvalue"=as.numeric(result$perm.pval[,ct]),
                                 "limma.category"=rep("TN", nrow(obs$obs.stat)),
                                 "perm.category"=rep("TN", nrow(obs$obs.stat))))
      curr_df$Amean<-as.numeric(curr_df$Amean)
      curr_df$coef<-as.numeric(curr_df$coef)
      curr_df$limma.category<-factor(curr_df$limma.category,
                                     levels = c("TN","TP","FP", "FN"))
      curr_df$limma.category[match(intersect(markers[[curr.ct.name]],
                                             row.names(result$obs.pval)[result$obs.pval.adj[,ct]<0.05]),  
                                   row.names(result$obs.pval))]<-"TP"
      curr_df$limma.category[match(setdiff(row.names(result$obs.pval)[result$obs.pval.adj[,ct]< 0.05], 
                                           markers[[curr.ct.name]]),
                                   row.names(result$obs.pval))]<-"FP"
      curr_df$limma.category[match(intersect(markers[[curr.ct.name]], 
                                             row.names(result$obs.pval)[result$obs.pval.adj[,ct]> 0.05]),
                                   row.names(result$obs.pval))]<-"FN"
      
      curr_df$perm.category<-factor(curr_df$perm.category, levels = c("TN","TP","FP", "FN"))
      curr_df$perm.category[match(intersect(markers[[curr.ct.name]], 
                                            row.names(result$obs.pval)[result$perm.pval.adj[,ct]<0.05]),  
                                  row.names(result$obs.pval))]<-"TP"
      curr_df$perm.category[match(setdiff(row.names(result$obs.pval)[result$perm.pval.adj[,ct]< 0.05], 
                                          markers[[curr.ct.name]]),
                                  row.names(result$obs.pval))]<-"FP"
      curr_df$perm.category[match(intersect(markers[[curr.ct.name]], 
                                            row.names(result$obs.pval)[result$perm.pval.adj[,ct]> 0.05]),
                                  row.names(result$obs.pval))]<-"FN"
      curr_df$limma.pvalue<-as.numeric(curr_df$limma.pvalue)
      curr_df$perm.pvalue<-as.numeric(curr_df$perm.pvalue)
      
      
      rb.df<- data.frame(genes=rep(row.names(curr_df), times=2),
                         pvalue =c(curr_df$limma.pvalue, curr_df$perm.pvalue),
                         category=c(curr_df$limma.category, curr_df$perm.category),
                         test.name=rep(c("limma","permutation"), each=nrow(curr_df))
      )
      rb.df$category<- factor(rb.df$category)
      rb.df$test.name<- factor(rb.df$test.name)
      rb.df$category<- factor(rb.df$category, levels = c(1,2,3,4),
                              labels =c("TN", "TP", "FP", "FN") )
      # raw p-value plot
      plot(ggplot(data=rb.df,
                  aes(x=pvalue, color=category,fill = test.name, 
                      y = test.name,shape=category, group = test.name))+
             geom_point(position=position_jitterdodge(jitter.width=2, jitter.height=0), 
                        size=2, alpha=0.8)+
             geom_hline(yintercept =1/perm.size, col="black")+
             xlab(paste("raw pvalue"))+
             scale_shape_manual(values=c(20, 8, 15, 17))+
             scale_color_manual(values=c("black","steelblue","red","orange"))+
             theme(legend.position= "bottom")+
             ggtitle(paste("raw pvalues for", curr.ct.name, "cells"))+
             ylab("")+guides(fill = "none")+
             facet_zoom(xlim=c(0, 0.05),horizontal=FALSE,split=F)
      )
      curr_df$true.marker.genes<- FALSE
      curr_df[match(intersect(markers[[curr.ct.name]], row.names(curr_df)), 
                    row.names(curr_df)), "true.marker.genes"]=TRUE
      curr_df$true.marker.genes<- factor(curr_df$true.marker.genes)
      table(curr_df$true.marker.genes)
      
      curr_df$true.marker.genes<- FALSE
      curr_df[match(intersect(markers[[curr.ct.name]], row.names(curr_df)), row.names(curr_df)), "true.marker.genes"]=TRUE
      curr_df$true.marker.genes<- factor(curr_df$true.marker.genes)
      table(curr_df$true.marker.genes)
      plot(ggplot(data=curr_df, aes(x=limma.pvalue, y=perm.pvalue,
                                    color=true.marker.genes, 
                                    size=true.marker.genes))+
             geom_point(
               #position=position_jitterdodge(jitter.width=0, 
               #                              jitter.height=0.001), 
               alpha=0.7)+
             ylab("permutation raw pvalue")+
             theme(legend.position= "bottom")+
             xlab("limma raw pvalue")+
             scale_color_manual(values = c("black", "steelblue"))+
             ggtitle(paste(curr.ct.name, "cells"))+
             geom_hline(yintercept = 0.05, color = "black")+
             geom_vline(xintercept = 0.05, color = "black")+
             geom_abline(slope=1, intercept = 0,lwd=1, color = "orange", linetype="dashed")+
             facet_zoom(xlim=c(0, 0.1),ylim=c(0, 0.1), horizontal = FALSE)+
             scale_size_manual(values = c(1, 2), drop=F))
      
      
      limma.pval <- sort(obs$obs.pval[,ct])
      expected.limma.pval <- 1:length(limma.pval)/length(limma.pval)
      
      perm.pval <- sort(result$perm.pval[,ct])
      expected.perm.pval <- 1:length(perm.pval)/length(perm.pval)
      
      pv_df<- data.frame(expected = expected.perm.pval, pval = perm.pval, 
                         true.DE=FALSE)
      pv_df[match(intersect(markers[[curr.ct.name]], row.names(pv_df)), 
                  row.names(pv_df)), "true.DE"]=TRUE
      plot(ggplot(pv_df, aes(x=-log10(expected), y=-log10(pval), color=true.DE)) + 
             geom_point( size = 1, alpha=0.7) + 
             scale_color_manual(values = c("black","steelblue"))+
             geom_line(aes(x=-log10(expected), y=-log10(expected)),
                       linetype="dashed", color = "black") +
             labs(y=expression(Permutation~~Observed~~-log[10](italic(p))), 
                  x=expression(Expected~~-log[10](italic(p))),
                  color="true marker genes") +
             ggtitle(paste(curr.ct.name,"cells",sep=" "))+
             theme(legend.position="bottom")
      )
      
      pv_df<- data.frame(expected = expected.limma.pval, pval = limma.pval, 
                         true.DE=FALSE)
      pv_df[match(intersect(markers[[curr.ct.name]], row.names(pv_df)), 
                  row.names(pv_df)), "true.DE"]=TRUE
      plot(ggplot(pv_df, aes(x=-log10(expected), y=-log10(pval), color=true.DE)) + 
             geom_point( size = 1, alpha=0.7) + 
             scale_color_manual(values = c("black","steelblue"))+
             geom_line(aes(x=-log10(expected), y=-log10(expected)),
                       linetype="dashed", color = "black") +
             labs(y=expression(Limma~~Observed~~-log[10](italic(p))), 
                  x=expression(Expected~~-log[10](italic(p))),
                  color="true marker genes") +
             ggtitle(paste(curr.ct.name,"cells",sep=" "))+
             theme(legend.position="bottom")
      )
      
    }else{
      perm.pval <- sort(result$perm.pval[,ct])
      expected.perm.pval <- 1:length(perm.pval)/length(perm.pval)
      
      pv_df<- data.frame(expected = expected.perm.pval, pval = perm.pval)
      plot(ggplot(pv_df, aes(x=-log10(expected), y=-log10(pval))) + 
             geom_point( size = 1, alpha=0.7) + 
             scale_color_manual(values = c("black","steelblue"))+
             geom_line(aes(x=-log10(expected), y=-log10(expected)),
                       linetype="dashed", color = "black") +
             ggtitle(paste(curr.ct.name,"cells",sep=" "))+
             labs(y=expression(Permutation~~Observed~~-log[10](italic(p))), 
                  x=expression(Expected~~-log[10](italic(p)))) 
      )
      
      limma.pval <- sort(obs$obs.pval[,ct])
      expected.limma.pval <- 1:length(limma.pval)/length(limma.pval)
      pv_df<- data.frame(expected = expected.limma.pval, pval = limma.pval)
      plot(ggplot(pv_df, aes(x=-log10(expected), y=-log10(pval))) + 
             geom_point( size = 1, alpha=0.7) + 
             scale_color_manual(values = c("black","steelblue"))+
             geom_line(aes(x=-log10(expected), y=-log10(expected)),
                       linetype="dashed", color = "black") +
             ggtitle(paste(curr.ct.name,"cells",sep=" "))+
             labs(y=expression(Limma~~Observed~~-log[10](italic(p))), 
                  x=expression(Expected~~-log[10](italic(p)))) 
      )
      
    }
    
    
  }
  dev.off()
}