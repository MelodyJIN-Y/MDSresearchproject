computeSummary<- function(obs.p, perm.p, true.upreg, cell_num){
  
  summ<- as.data.frame(matrix(NA, nrow =8 , ncol =ncol(obs.p)))
  colnames(summ)<- colnames(obs.p)
  row.names(summ)<- c("#cells", "#upreg(limma)", "#trueUpreg(limma)", "%trueUpreg(limma)",
                         "#upreg(permutation)","#trueUpreg(permutation)", "%trueUpreg(permutation)",
                         "#overlap(limma&permutation)")
  for (i in 1: ncol(obs.p)){
    limma_upreg_genes =  row.names(obs.p[which(obs.p[,i]<0.05),])
    perm_upreg_genes =  row.names(perm.p[which(perm.p[,i]<0.05),])
    summ[1,i]<- cell_num[i]
    summ[2,i]<-length(limma_upreg_genes)
    summ[3,i]<-length(intersect(true.upreg[[i]], limma_upreg_genes))
    summ[4,i]<-100*round(length(intersect(true.upreg[[i]], limma_upreg_genes)) / length(true.upreg[[i]]),5)
    summ[5,i]<-length(perm_upreg_genes)
    summ[6,i]<-length(intersect(true.upreg[[i]], perm_upreg_genes))
    summ[7,i]<-100*round(length(intersect(true.upreg[[i]], perm_upreg_genes)) / length(true.upreg[[i]]),5)
    summ[8,i]<-length(intersect(limma_upreg_genes, perm_upreg_genes))
    
  }
  
  return (summ)
}