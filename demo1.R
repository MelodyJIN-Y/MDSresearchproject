cm<- matrix(rbinom(n=24, prob = 0.2, size = 4),nrow = 3, ncol = 8)
cm
obs.stat<- t(apply(cm, 1, function(r) c(t.test(r[c(1,2)],r[-c(1,2)])$statistic, 
             t.test(r[c(3,4,5)],r[-c(3,4,5)])$statistic,
             t.test(r[c(6,7,8)],r[-c(6,7,8)])$statistic)))
#set.seed(10)
perm.stat<- matrix(NA, nrow = 100, ncol = 3)

for (i in 1:100) {
  shuffle <- sample(8, replace = FALSE, size =8)
  cm.perm<- cm[, shuffle]
  perm.stat[i, ]<- t(apply(cm.perm, 1, function(r) c(t.test(r[c(1,2)],r[-c(1,2)])$statistic, 
                                                t.test(r[c(3,4,5)],r[-c(3,4,5)])$statistic,
                                                t.test(r[c(6,7,8)],r[-c(6,7,8)])$statistic)))

  
}
 