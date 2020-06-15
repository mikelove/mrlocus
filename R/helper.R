#' Collapse high correlation SNPs
#'
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#' @param threshold threshold on absolute value of
#' correlation for collapsing, e.g. will collapse
#' if SNPs are more correlated (or anti-correlated)
#' than this amount
#' @param plot logical, draw a before/after grid of plots
#'
#' @return list with modified ld_mat and sum_stat lists
#' 
#' @export
collapseHighCorSNPs <- function(sum_stat, ld_mat, thresh=.95, plot=TRUE) {
  stopifnot(length(sum_stat) == length(ld_mat))
  nsnps <- formatC(sapply(sum_stat,nrow), width=2, flag=0)
  message(paste0("pre:  ",paste(nsnps,collapse=",")))
  gs <- list()
  for (j in seq_along(sum_stat)) {
    if (nrow(sum_stat[[j]]) == 1) next
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
  nsnps <- formatC(sapply(sum_stat,nrow), width=2, flag=0)
  message(paste0("post: ",paste(nsnps,collapse=",")))
  list(sum_stat=sum_stat, ld_mat=ld_mat)
}

#' Flip alleles and gather results into lists
#'
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#' @param a name of A in columns of sum_stat ("eQTL")
#' @param b name of B ("GWAS")
#' @param ref name of reference allele
#' @param eff name of effect allele
#' @param beta name of estimated coefficient
#' @param se name of standard error
#' @param major_plink name of major plink allele
#' @param sep character separator in column names
#' @param ab_last A/B descriptor is last in column names
#' (e.g. "beta_eqtl", "se_eqtl"))
#' @param plot logical, draw a scatterplot of the flipped betas
#'
#' @return list with estimated coefficients, standard
#' errors, LD matrix, and allele table
#' 
#' @export
flipAllelesAndGather <- function(sum_stat, ld_mat,
                                 a, b, ref, eff,
                                 beta, se, major_plink,
                                 sep, ab_last=TRUE, plot=TRUE) {

  # the following allow for arbitrary incoming column names.
  # the point of this is to reduce mistakes that might occur
  # if users manually had to modify their column names.
  ref_a_nm <- if (ab_last) paste(ref, a, sep=sep) else paste(a, ref, sep=sep)
  ref_b_nm <- if (ab_last) paste(ref, b, sep=sep) else paste(b, ref, sep=sep)
  eff_a_nm <- if (ab_last) paste(eff, a, sep=sep) else paste(a, eff, sep=sep)
  eff_b_nm <- if (ab_last) paste(eff, b, sep=sep) else paste(b, eff, sep=sep)
  beta_a_nm <- if (ab_last) paste(beta, a, sep=sep) else paste(a, beta, sep=sep)
  beta_b_nm <- if (ab_last) paste(beta, b, sep=sep) else paste(b, beta, sep=sep)
  se_b_nm <- if (ab_last) paste(se, b, sep=sep) else paste(b, se, sep=sep)
  se_a_nm <- if (ab_last) paste(se, a, sep=sep) else paste(a, se, sep=sep)
  
  beta_hat_a <- list()
  beta_hat_b <- list()
  se_a <- list()
  se_b <- list()
  Sigma <- list()
  alleles <- list()
  
  for (j in seq_along(sum_stat)) {
    # the index SNP in the signal cluster for A (eQTL) based on Z stat
    idx <- which.max(abs(sum_stat[[j]][[beta_a_nm]]/sum_stat[[j]][[se_a_nm]]))
    # check: reference B (GWAS) allele must be either reference or effect allele in A (eQTL)
    stopifnot(sum_stat[[j]][[ref_b_nm]] == sum_stat[[j]][[ref_a_nm]] |
              sum_stat[[j]][[ref_b_nm]] == sum_stat[[j]][[eff_a_nm]])
    # flip B (GWAS) so that same ref allele is described for B (GWAS) as for A (eQTL)
    flip <- which(sum_stat[[j]][[ref_b_nm]] != sum_stat[[j]][[ref_a_nm]])
    sum_stat[[j]]$beta_b_flipped <- sum_stat[[j]][[beta_b_nm]]
    sum_stat[[j]]$beta_b_flipped[flip] <- -1 * sum_stat[[j]][[beta_b_nm]][flip]
    # flip alleles other than the index so they have positive correlation (LD) with index
    # and based on what the major allele is according to plink
    ld.sign <- sign(ld_mat[[j]][,idx])
    plink.agree <- ifelse(sum_stat[[j]][[major_plink]] == sum_stat[[j]][[ref_a_nm]], 1, -1)
    beta_a <- plink.agree * ld.sign * sum_stat[[j]][[beta_a_nm]]
    beta_b <- plink.agree * ld.sign * sum_stat[[j]]$beta_b_flipped
    # only flip LD matrix based on positive correlation with index (bc it comes from plink)
    ld.flipped <- t(t(ld_mat[[j]]) * ld.sign) * ld.sign
    Sigma[[j]] <- ld.flipped
    # record the alleles after all the flipping
    alleles[[j]] <- data.frame(ref=sum_stat[[j]][[ref_a_nm]],
                               eff=sum_stat[[j]][[eff_a_nm]])
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
    se_a[[j]] <- sum_stat[[j]][[se_a_nm]]
    se_b[[j]] <- sum_stat[[j]][[se_b_nm]]
  }
  out <- list(beta_hat_a=beta_hat_a,
              beta_hat_b=beta_hat_b,
              se_a=se_a, se_b=se_b,
              Sigma=Sigma,
              alleles=alleles)
  if (plot) {
    plotInitEstimates(out, a=a, b=b)
  }
  return(out)
}

#' Plot initial estimates over signal clusters
#'
#' @param x list of signal clusters data with beta_hat_a
#' and beta_hat_b lists
#' @param a name of A experiment
#' @param b name of B experiment
#'
#' @export
plotInitEstimates <- function(x, a="eQTL", b="GWAS") {
  nsnp <- lengths(x$beta_hat_a)
  plot(unlist(x$beta_hat_a), unlist(x$beta_hat_b),
       xlab=paste("beta", a),
       ylab=paste("beta", b), 
       col=rep(seq_along(nsnp),nsnp),
       pch=rep(seq_along(nsnp),nsnp))
  text(unlist(x$beta_hat_a), unlist(x$beta_hat_b),
       do.call(c, sapply(nsnp, seq_len)), pos=4, cex=.5,
       col=rep(seq_along(nsnp),nsnp))
  abline(h=0, col=rgb(0,0,0,.3))
}
