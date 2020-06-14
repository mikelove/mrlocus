#' Collapse high correlation SNPs
#'
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#' @param threshold threshold on absolute value of
#' correlation for collapsing, e.g. will collapse
#' if SNPs are more correlated (or anti-correlated)
#' than this amount
#' @param plot to draw a before after grid of plots
#'
#' @return list with modified ld_mat and sum_stat lists
#' 
#' @export
collapseHighCorSNPs <- function(sum_stat, ld_mat, thresh=.95, plot=TRUE) {
  stopitnot(length(sum_stat) == length(ld_mat))
  gs <- list()
  for (j in seq_along(sum_stat)) {
    gs[[2*j-1]] <- pheatmap::pheatmap(ld_mat[[j]], breaks=seq(-1,1,length=101),
                                      cluster_rows=FALSE, cluster_cols=FALSE,
                                      show_rownames=FALSE, show_colnames=FALSE,
                                      border_color=NA, 
                                      silent=TRUE)$gtable
    hc <- hclust(as.dist(1-abs(ld_mat[[j]])))
    #plot(hc); abline(h=1-thresh,col="red")
    hclusters <- cutree(hc, h=1-thresh)
    nhclust <- max(hclusters)
    # greedy reduction: just take the first SNP
    # for each hierarchical cluster
    idx <- match(seq_len(nhclust), hclusters) 
    ld_mat[[j]] <- ld_mat[[j]][idx,idx]
    sum_stat[[j]] <- sum_stat[[j]][idx,]
    gs[[2*j]] <- pheatmap::pheatmap(ld_mat[[j]], breaks=seq(-1,1,length=101),
                                    cluster_rows=FALSE, cluster_cols=FALSE,
                                    show_rownames=FALSE, show_colnames=FALSE,
                                    border_color=NA, 
                                    silent=TRUE)$gtable
  }
  if (plot) {
    gridExtra::grid.arrange(ncol=2, grobs=gs)
  }
  list(sum_stat=sum_stat, ld_mat=ld_mat)
}

#' Flip alleles and gather results into lists
#'
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#'
#' @return list with estimated coefficients, standard
#' errors, LD matrix, and allele table
#' 
#' @export
flipAllelesAndGather <- function(sum_stat, ld_mat) {
  beta_hat_a <- list()
  beta_hat_b <- list()
  se_a <- list()
  se_b <- list()
  Sigma <- list()
  alleles <- list()
  for (j in seq_along(sum_stat)) {
    # the index SNP in the signal cluster for A (eQTL) based on Z stat
    idx <- which.max(abs(sum_stat[[j]]$beta_eQTL/sum_stat[[j]]$se_eQTL))
    # check: reference B (GWAS) allele must be either reference or effect allele in A (eQTL)
    stopifnot(sum_stat[[j]]$Ref_GWAS == sum_stat[[j]]$Ref_eQTL |
              sum_stat[[j]]$Ref_GWAS == sum_stat[[j]]$Effect_eQTL)
    # flip B (GWAS) so that same ref allele is described for B (GWAS) as for A (eQTL)
    flip <- which(sum_stat[[j]]$Ref_GWAS != sum_stat[[j]]$Ref_eQTL)
    sum_stat[[j]]$beta_GWAS_flipped <- sum_stat[[j]]$beta_GWAS
    sum_stat[[j]]$beta_GWAS_flipped[flip] <- -1 * sum_stat[[j]]$beta_GWAS[flip]
    # flip alleles other than the index so they have positive correlation (LD) with index
    # and based on what the major allele is according to plink
    ld.sign <- sign(ld_mat[[j]][,idx])
    plink.agree <- ifelse(sum_stat[[j]]$Major_plink == sum_stat[[j]]$Ref_eQTL, 1, -1)
    beta_a <- plink.agree * ld.sign * sum_stat[[j]]$beta_eQTL
    beta_b <- plink.agree * ld.sign * sum_stat[[j]]$beta_GWAS_flipped
    # only flip LD matrix based on positive correlation with index (bc it comes from plink)
    ld.flipped <- t(t(ld_mat[[j]]) * ld.sign) * ld.sign
    Sigma[[j]] <- ld.flipped
    # record the alleles after all the flipping
    alleles[[j]] <- data.frame(ref=sum_stat[[j]]$Ref_eQTL,
                               eff=sum_stat[[j]]$Effect_eQTL)
    idx2 <- ld.sign * plink.agree == -1
    tmp.ref <- alleles[[j]]$ref[idx2]
    tmp.eff <- alleles[[j]]$eff[idx2]
    alleles[[j]]$ref[idx2] <- tmp.eff
    alleles[[j]]$eff[idx2] <- tmp.ref
    # finally, flip alleles so that effect size is positive for index SNP for A (eQTL)
    if (beta_a[idx] < 0) {
      beta_a <- beta_a * -1
      beta_b <- beta_b * -1
      tmp.ref <- alleles[[j]]$ref
      tmp.eff <- alleles[[j]]$eff
      alleles[[j]]$ref <- tmp.eff
      alleles[[j]]$eff <- tmp.ref
    }
    beta_hat_a[[j]] <- beta_a
    beta_hat_b[[j]] <- beta_b
    se_a[[j]] <- sum_stat[[j]]$se_eQTL
    se_b[[j]] <- sum_stat[[j]]$se_GWAS
  }
  return(list(beta_hat_a=beta_hat_a,
              beta_hat_b=beta_hat_b,
              se_a=se_a, se_b=se_b,
              Sigma=Sigma,
              alleles=alleles))
}
