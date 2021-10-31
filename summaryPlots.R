# This is a helper function which will generate MA and ROC plots
# Input:  plt.title = file name, a string
#         result = An object returned by function runPermuation(), a list
#         obs = An object returned by function computeObservation() ,  a list
#         statistic = used statistic, a string
#         hlins = a number used to add a horizontal line to the MA plot

summaryPlots<- function(plt.title, result, markers, obs, statistic, hlins = 0){
  pdf(paste(plt.title, "pdf", sep="."))
  auc.values<- as.data.frame(matrix(0, nrow=2, ncol=ncol(result$obs.pval)+1))
  colnames(auc.values)<- c("test", colnames(result$obs.pval))
  auc.values$test<- c("limma","permutation")
  for (ct in 1:ncol(result$obs.pval)){
    curr.ct.name<- colnames(result$obs.pval)[ct]
    if (curr.ct.name %in% names(markers)){
      
      # MA plot
      curr_df<- data.frame(cbind("Amean"=obs$fit.cont.eb.obs$Amean, 
                                 "coef"=obs$fit.cont.eb.obs$coefficients[,ct],
                                 "limma.pvalue"=as.numeric(result$obs.pval.adj[,ct]),
                                 "perm.pvalue"=as.numeric(result$perm.pval.adj[,ct]),
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
     # breaks<-c(0.1,0.2, 0.3, 0.4, 0.5, 0.7, 0.95, 0.99)
    #  labels<- c("0.9","0.8","0.7","0.6", "0.5", "0.3", "0.05", "0.01" )
      limma.p<- ggplot(data=curr_df, 
                       aes(x =Amean, y=coef, color = limma.category, 
                           shape= limma.category, size=limma.category
                       ))+
        geom_point(position=position_jitterdodge(jitter.width=0, jitter.height=0.002), alpha=0.8)+
        geom_hline(yintercept=hlins, size=1, color = "darkgreen")+
        geom_hline(yintercept=0, size=1, color = "black" )+
        scale_color_manual(values = c("black","royalblue2","red","goldenrod2"),drop=F)+
        geom_text(aes(-1.5,hlins,label = paste(hlins),hjust=1,  vjust = -0.5), size=3, color="darkgreen")+
        theme(
          legend.title=element_blank(),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15), 
          panel.spacing = unit(0.5, "lines"), 
          legend.position="bottom",
          legend.text=element_text(margin=margin(r=0.3,unit="cm"), size=15))+
        xlab("Average log-expression")+ylab("log-fold-change")+
        ggtitle(paste("limma", statistic))+
        scale_y_continuous(breaks=c(-1, -2, 0, 1, 2 ))+
        #scale_color_gradientn(name="Legend 2",
        #                     #labels=comma, 
        #                     limits=c(0, 1),
        #                     #values = c(0, 0.5, 1),
        #                     colours=rev(c("#000000", "#FFFFFF", "#BA0000")), 
        #                     values=c(0, 0.053, 1)) +
        scale_shape_manual(values=c(20, 8, 15, 17),drop=F)+
        scale_size_manual(values=c(0.1, 1.5,1.5,1.5),drop=F)
      
      perm.p<- ggplot(data=curr_df, 
                      aes(x =Amean, y=coef, color = perm.category, 
                          shape= perm.category,size=perm.category
                          ))+
        geom_point(position=position_jitterdodge(jitter.width=0, jitter.height=0.002), alpha=0.9)+
        geom_hline(yintercept=hlins, size=1, color = "darkgreen")+
        geom_hline(yintercept=0, size=1, color = "black" )+
        geom_text(aes(-1.5,hlins,label = paste(hlins),hjust=1,  vjust = -0.5), size=3, color="darkgreen")+
        scale_color_manual(values = c("black","royalblue2","red","goldenrod2"),drop=F)+
        ggtitle(paste("permutation", statistic))+
        theme(
          legend.title=element_blank(),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15), 
          panel.spacing = unit(0.5, "lines"), 
          legend.position="bottom",
          legend.text=element_text(margin=margin(r=0.3,unit="cm"), size=15))+
        xlab("Average log-expression")+ylab("log-fold-change")+
        scale_shape_manual(values=c(20, 8, 15, 17),drop=F)+
        scale_y_continuous(breaks=c(-1, -2, 0, 1, 2 ))+
        scale_size_manual(values=c(0.1, 1.5,1.5,1.5),drop=F)
      pl<- ggarrange(limma.p, perm.p,  ncol =1,  nrow = 2,
                     legend="bottom",align="v", common.legend=T)
      plot(annotate_figure(pl, 
                      top = text_grob(paste(curr.ct.name, "cells"), 
                                      color = "black", face = "bold", size = 14)))
      
      # AUC plot
      resp=list()
      
        resp.g <- as.vector(ifelse(test=(row.names(result$perm.pval.adj) %in% markers[[curr.ct.name]]),  
                                   yes=1, no=0))
        resp[[curr.ct.name]]<-resp.g

      

        
        pred.perm<-result$perm.pval.adj[,curr.ct.name]
        pred.limma<-result$obs.pval.adj[,curr.ct.name]
        
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
        
        
        rb.df<- data.frame(genes=rep(row.names(curr_df), times=2),
          pvalue =c(curr_df$limma.pvalue, curr_df$perm.pvalue),
                              category=c(curr_df$limma.category, curr_df$perm.category),
                           test.name=rep(c("limma","permutation"), each=nrow(curr_df))
                              )
        rb.df$category<- factor(rb.df$category)
        rb.df$test.name<- factor(rb.df$test.name)
        rb.df$category<- factor(rb.df$category, levels = c(1,2,3,4),
                                labels =c("TN", "TP", "FP", "FN") )
        # adjusted p-value plot
        plot(ggplot(data=rb.df,
                          aes(x=pvalue, color=category,fill = test.name, 
                              y = test.name,shape=category, group = test.name))+
               geom_point(position=position_jitterdodge(jitter.width=2, jitter.height=0), 
                          size=2, alpha=0.8)+
               geom_hline(yintercept =1/perm.size, col="black")+
               xlab(paste("adjusted pvalue"))+
               scale_shape_manual(values=c(20, 8, 15, 17))+
               scale_color_manual(values=c("black","steelblue","red","orange"))+
               theme(legend.position= "bottom")+
               ggtitle(paste("adjusted pvalues for", curr.ct.name, "cells"))+
               ylab("")+guides(fill = "none")+
               facet_zoom(xlim=c(0, 0.05),horizontal=FALSE,split=F)
        )
        curr_df$true.marker.genes<- FALSE
        curr_df[match(intersect(markers[[curr.ct.name]], row.names(curr_df)), row.names(curr_df)), "true.marker.genes"]=TRUE
        curr_df$true.marker.genes<- factor(curr_df$true.marker.genes)
        table(curr_df$true.marker.genes)
        plot(ggplot(data=curr_df, aes(x=limma.pvalue, y=perm.pvalue,
                                      color=true.marker.genes, 
                                      size=true.marker.genes))+
               geom_point(
                           alpha=0.7)+
               ylab("permutation adjusted pvalue")+
               theme(legend.position= "bottom")+
               xlab("limma adjusted pvalue")+
               scale_color_manual(values = c("black", "steelblue"))+
               ggtitle(paste(curr.ct.name, "cells"))+
               geom_hline(yintercept = 0.05, color = "black")+
               geom_vline(xintercept = 0.05, color = "black")+
               geom_abline(slope=1, intercept = 0,lwd=1, color = "orange", linetype="dashed")+
               facet_zoom(xlim=c(0, 0.1),ylim=c(0, 0.1), horizontal = FALSE)+
          scale_size_manual(values = c(1, 2), drop=F))
        
    
       
    }else{
      
      plot(obs$fit.cont.eb.obs$Amean, obs$fit.cont.eb.obs$coefficients[,ct], 
           pch=16, cex=0.8, main = paste(curr.ct.name,"cells (limma)",sep=" "), 
           xlab="Average log-expression", ylab="log-fold-change")
      points(obs$fit.cont.eb.obs$Amean[result$obs.pval.adj[,ct]<0.05], 
             obs$fit.cont.eb.obs$coefficients[result$obs.pval.adj[,ct]<0.05,ct], 
             pch=16, col="steelblue")
      abline(h=0)
      legend(x = "bottomright",          
             legend = c("DEs", "non-DEs"),
             col = c("steelblue", "black"),           
             lwd = 2)  
      
      plot(obs$fit.cont.eb.obs$Amean, obs$fit.cont.eb.obs$coefficients[,ct], 
           pch=16, cex=0.8, main = paste(curr.ct.name,"cells (permutation)",sep=" "), 
           xlab="Average log-expression", ylab="log-fold-change")
      points(obs$fit.cont.eb.obs$Amean[result$perm.pval.adj[,ct]<0.05], 
             obs$fit.cont.eb.obs$coefficients[result$perm.pval.adj[,ct]<0.05,ct], 
             pch=16, col="steelblue")
      abline(h=0)
      legend(x = "bottomright",          
             legend = c("DEs", "non-DEs"),
             col = c("steelblue", "black"),           
             lwd = 2) 

      
    }
    
    
  }
  dev.off()
  return (auc.values)
}