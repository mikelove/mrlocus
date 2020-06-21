context("collapse")
library(mrlocus)

test_that("collapsing high cor SNPs", {

  sum_stat <- list(data.frame(beta.eqtl=c(1,2,.5,1,2),
                              se.eqtl=c(.5,.5,.5,.5,.5),
                              snp=letters[1:5]))
  sum_stat <- lapply(sum_stat, function(x) {
    x$abs.z <- abs(x$beta.eqtl / x$se.eqtl)
    x
  })

  ld <- matrix(c(1,.99,0,0,0,
                 .99,1,0,0,0,
                 0,0,1,0,0,
                 0,0,0,1,.99,
                 0,0,0,.99,1
                 ), byrow=TRUE, ncol=5)
  ld_mat <- list(ld)

  out <- collapseHighCorSNPs(sum_stat, ld_mat)
  out$sum_stat[[1]]
  
  out <- collapseHighCorSNPs(sum_stat, ld_mat, score="abs.z")
  out$sum_stat[[1]]
  
})
