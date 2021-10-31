computeAUC<- function(result, markers){
  auc.values<- as.data.frame(matrix(0, nrow=2, ncol=ncol(result$obs.pval)+1))
  colnames(auc.values)<- c("test", colnames(result$obs.pval))
  auc.values$test<- c("limma","permutation")
  for (ct in 1:ncol(result$obs.pval)){
    curr.ct.name<- colnames(result$obs.pval)[ct]
    if (curr.ct.name %in% names(markers)){
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
    }else{ 
      auc.values[auc.values$test=="limma",curr.ct.name]<- 0
      auc.values[auc.values$test=="permutation",curr.ct.name]<- 0
    }
  }
  return (auc.values)
}
