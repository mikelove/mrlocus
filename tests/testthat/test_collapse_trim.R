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
  ld_mat2 <- ld_mat

  sum_stat[[1]]
  
  out <- collapseHighCorSNPs(sum_stat, ld_mat, ld_mat2)
  dev.off()
  out$sum_stat[[1]]
  
  out <- collapseHighCorSNPs(sum_stat, ld_mat, ld_mat2, score="abs.z")
  dev.off()
  out$sum_stat[[1]]
  
})


test_that("trimming clusters works", {


  r2 <- matrix(0, ncol=3, nrow=3)
  r2[1,2] <- r2[2,1] <- .2
  diag(r2) <- 1
  trim_clusters <- trimClusters(r2, r2_threshold=0.05) # should tell us to trim #2
  expect_equal(trim_clusters, 2)

})
