# A function to geenrate plots to compare Limma and permutation performance 
# marker.genes is a list of evident marker genes for each cell type, in the same order as the columns of obs.p
# if marker.genes is NULL, 
# the columns "TP", "precision","sensitivity","F1","TP.overlap","AUC" will all be zero

# unique.marker.genes is a prameter to specify the length of unique markers for each cell type
# this is included to address the gene synnomym issues read from the marker gene dataframe 
# if unique.marker.genes = NULL (the defualt), all genes in marker genes will be assumed to be unique for each cell type 

computeComparison<- function(obs.p, perm.p, marker.genes, 
                             unique.marker.genes=NULL, ids,
                             cell_num, pvalue.cutoff=0.05, 
                             perm.size, auc.values,statistic){

    summ<- as.data.frame(matrix(NA, nrow =ncol(obs.p)*2 , ncol =17))
    colnames(summ)<- c("id","cell.type","num.cells", "total.num.cells","test.name", 
                       "num.true.marker","num.upreg","TP", "precision", 
                       "sensitivity","F1", "overlap","TP.overlap","pvalue.cutoff",
                       "perm.size", "AUC","statistic")
    summ[, "cell.type"]<- rep(colnames(obs.p), 2)
    summ[, "id"]<- ids
    summ[, "test.name"]<- "permutation"
    summ[1:ncol(obs.p), "test.name"]<- rep("limma", ncol(obs.p))
    summ[, "statistic"]<- statistic
    summ[, "perm.size"]<- perm.size
    summ[, "total.num.cells"]<-sum(cell_num)
    summ[, "pvalue.cutoff"]<- pvalue.cutoff
    
    
    for (i in 1: ncol(obs.p)){
      summ[i, "num.cells"]<- cell_num[i]
      summ[i+ncol(obs.p), "num.cells"]<- cell_num[i]
      
      summ[i, "num.true.marker"]<- length(unique.marker.genes[[colnames(obs.p)[i]]])
      summ[i+ncol(obs.p), "num.true.marker"]<- length(unique.marker.genes[[colnames(obs.p)[i]]])
      
      # positive genes 
      limma_upreg_genes =  row.names(obs.p[which(obs.p[,i] < pvalue.cutoff),])
      perm_upreg_genes =  row.names(perm.p[which(perm.p[,i] < pvalue.cutoff),])
      
      # negative genes 
      limma_nonsig_genes =  setdiff(row.names(obs.p), limma_upreg_genes)
      perm_nonsig_genes =  setdiff(row.names(obs.p), perm_upreg_genes)
      
      overlap<-length(intersect(limma_upreg_genes, perm_upreg_genes))
      summ[i, "num.upreg"]<- length(limma_upreg_genes)
      summ[i+ncol(obs.p), "num.upreg"]<- length(perm_upreg_genes)
      summ[i, "overlap"]<-overlap
      summ[i+ncol(obs.p), "overlap"]<- overlap
      
      if (is.null(marker.genes[[colnames(obs.p)[i]]]) ==TRUE){
        
        summ[i, c("TP", "precision","sensitivity","F1","TP.overlap","AUC")]<- 0
        summ[i+ncol(obs.p), c("TP", "precision","sensitivity","F1","TP.overlap","AUC")]<- 0
        
      }else{
        # metrics 
        tp.l<- length(intersect(limma_upreg_genes, marker.genes[[colnames(obs.p)[i]]]))
        fp.l<- length(limma_upreg_genes)-length(intersect(marker.genes[[colnames(obs.p)[i]]], limma_upreg_genes))
        tn.l<- length(intersect(setdiff(row.names(obs.p),marker.genes[[colnames(obs.p)[i]]]),limma_nonsig_genes))
        fn.l<- length(setdiff(marker.genes[[colnames(obs.p)[i]]], limma_upreg_genes))
        
        tp.p<- length(intersect(perm_upreg_genes, marker.genes[[colnames(obs.p)[i]]]))
        fp.p<- length(perm_upreg_genes)-length(intersect(marker.genes[[colnames(obs.p)[i]]], perm_upreg_genes))
        tn.p<- length(intersect(setdiff(row.names(obs.p),marker.genes[[colnames(obs.p)[i]]]),perm_nonsig_genes))
        fn.p<- length(setdiff(marker.genes[[colnames(obs.p)[i]]], perm_upreg_genes))
        
        tp.overlap<-length(intersect(intersect(limma_upreg_genes, marker.genes[[colnames(obs.p)[i]]]), 
                                     intersect(perm_upreg_genes, marker.genes[[colnames(obs.p)[i]]])))
        
        precision.p<- round(tp.p/length(perm_upreg_genes),5)
        precision.l<- round(tp.l/length(limma_upreg_genes),5)
        # recall=TruePositives / (TruePositives + FalseNegatives)
        if (is.null(unique.marker.genes) == TRUE){
          unique.marker.genes = marker.genes
        }
        sensitivity.p<-  round(tp.p/length(unique.marker.genes[[colnames(obs.p)[i]]]),5)
        sensitivity.l<-  round(tp.l/length(unique.marker.genes[[colnames(obs.p)[i]]]),5)
        summ[i, "TP"]<- tp.l
        summ[i+ncol(obs.p), "TP"]<- tp.p
        
        summ[i, "precision"]<-precision.l
        summ[i+ncol(obs.p), "precision"]<- precision.p
        
        summ[i, "sensitivity"]<- sensitivity.l
        summ[i+ncol(obs.p), "sensitivity"]<- sensitivity.p
        
        summ[i, "F1"]<- round((2 * precision.l * sensitivity.l) / (precision.l + sensitivity.l), 5)
        summ[i+ncol(obs.p), "F1"]<- round((2 * precision.p * sensitivity.p) / (precision.p + sensitivity.p), 5)
        
        
        summ[i, "TP.overlap"]<- tp.overlap
        summ[i+ncol(obs.p), "TP.overlap"]<- tp.overlap
        
        summ[i, "AUC"]<- auc.values[auc.values$test=="limma", colnames(obs.p)[i]]
        summ[i+ncol(obs.p), "AUC"]<- auc.values[auc.values$test=="permutation", colnames(obs.p)[i]]
      }
    }
  
  
  return (summ)
}
