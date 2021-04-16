(cm<- matrix(rnorm(n=80,  mu = 10,sd=1 ),nrow = 10, ncol = 8))
cm
(obs.stat<- t(apply(cm, 1, 
                   function(r) c(t.test(r[c(1,2)],r[-c(1,2)], alternative = "greater")$statistic, 
                                 t.test(r[c(3,4,5)],r[-c(3,4,5)],alternative = "greater")$statistic,
                                 t.test(r[c(6,7,8)],r[-c(6,7,8)],alternative = "greater")$statistic))))

#set.seed(10)
perm.size<- 2000
(perm.array<- array(0, dim = c(10, 3, perm.size)))
for (i in 1:perm.size) {
  shuffle <- sample(8, replace = FALSE, size =8)
  cm.perm<- cm[, shuffle]
  perm.array[,,i]<- t(apply(cm.perm, 1, 
                            function(r) c(t.test(r[c(1,2)],r[-c(1,2)], alternative = "greater")$statistic, 
                                                t.test(r[c(3,4,5)],r[-c(3,4,5)],alternative = "greater")$statistic,
                                                t.test(r[c(6,7,8)],r[-c(6,7,8)],alternative = "greater")$statistic)))

}
# 

perm.res<-apply(expand.grid(x = 1:10, y = 1:3), 1, 
                function(r) (sum(perm.array[r[1],r[2], ] >obs.stat[r[1],r[2]])+1)/(perm.size+1) )
(perm.res<- matrix(perm.res, nrow=10, ncol = 3, byrow=FALSE))

(obs.p<- t(apply(cm, 1, 
                function(r) c(t.test(r[c(1,2)],r[-c(1,2)],alternative = "greater")$p.value, 
                              t.test(r[c(3,4,5)],r[-c(3,4,5)],alternative = "greater")$p.value,
                              t.test(r[c(6,7,8)],r[-c(6,7,8)],alternative = "greater")$p.value))))
