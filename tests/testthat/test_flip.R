context("flip")
library(mrlocus)

test_that("flip alleles", {

  set.seed(1)
  sum_stat <- list(data.frame(beta.eqtl=c(rnorm(3,1,.05),-2,-1),
                              se.eqtl=  c(1,1,1,1,1),
                              beta.gwas=c(rnorm(3,1,.05),-2,1),
                              se.gwas=  c(1,1,1,1,1),
                              ref.eqtl=c("A","A","A","A","A"),
                              eff.eqtl=c("T","T","T","T","T"),
                              ref.gwas=c("A","A","A","A","T"),
                              eff.gwas=c("T","T","T","T","A"),
                              snp=letters[1:5]))

  ld <- matrix(c(1,.5,.5,-.5,-.5,
                 .5,1,.5,-.5,-.5,
                 .5,.5,1,-.5,-.5,
                 -.5,-.5,-.5,1,.5,
                 -.5,-.5,-.5,.5,1
                 ), byrow=TRUE, ncol=5)
  ld_mat <- list(ld)
  ld_mat2 <- ld_mat

  out <- flipAllelesAndGather(sum_stat, ld_mat, ld_mat2,
                              a="eqtl", b="gwas",
                              ref="ref", eff="eff",
                              beta="beta", se="se",
                              a2_plink="ref.eqtl",
                              a2_plink_mat2="ref.eqtl",
                              snp_id="snp", sep=".")
  dev.off()

  names(out)

})
