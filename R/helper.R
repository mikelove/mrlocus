#' Collapse high correlation SNPs
#'
#' A helper function to collapse sets of highly correlated
#' SNPs within signal clusters. This is recommended to run
#' before \code{\link{flipAllelesAndGather}}, and before
#' \code{\link{fitBetaColoc}}.
#' 
#' @param sum_stat list of summary statistic tables,
#' which is a list over signal clusters.
#' Each element of the list should be a data.frame
#' describing the eQTL and GWAS summary statistics.
#' The only column in \code{sum_stat} that is used
#' by the function is \code{score} (optional)
#' @param ld_mat list of LD matrices across
#' signal clusters
#' @param ld_mat2 optional second list of LD matrices
#' (for different populations). it will be returned
#' alongside the first \code{ld_mat}, which is used
#' for the collapsing. The second list of LD
#' matrices is just subset to the same set of SNPs
#' as the first
#' @param threshold threshold on absolute value of
#' correlation for collapsing, e.g. will collapse
#' if SNPs are more correlated (or anti-correlated)
#' than this amount
#' @param score name of a column of sum_stat data.frames
#' with a score, such that collapsing will choose the
#' highest score SNP per collapsed cluster. Otherwise,
#' if set to NULL, the first SNP will be used
#' @param plot logical, draw a before/after grid of plots
#' @param snp_id name of SNP id column in \code{sum_stat},
#' if specified will output a column \code{collapsed}
#' that lists which SNP ids are represented in the output
#' (i.e. which highly correlated SNPs were collapsed).
#'
#' @return list with subset \code{ld_mat} and \code{sum_stat}
#' lists (and \code{ld_mat2} if provided)
#' 
#' @export
collapseHighCorSNPs <- function(sum_stat, ld_mat, ld_mat2=NULL,
                                threshold=.95, score=NULL, plot=TRUE,
                                snp_id=NULL) {
  stopifnot(length(sum_stat) == length(ld_mat))
  stopifnot(all(sapply(ld_mat, is, "matrix")))
  stopifnot(nrow(ld_mat[[1]]) == ncol(ld_mat[[1]]))
  stopifnot(nrow(ld_mat[[1]]) == nrow(sum_stat[[1]]))
  two.ld <- !is.null(ld_mat2)
  if (two.ld) {
    stopifnot(length(sum_stat) == length(ld_mat2))
    stopifnot(all(sapply(ld_mat2, is, "matrix")))
    stopifnot(nrow(ld_mat2[[1]]) == ncol(ld_mat2[[1]]))
    stopifnot(nrow(ld_mat2[[1]]) == nrow(sum_stat[[1]]))
  }
  if (plot) {
    if (!requireNamespace("pheatmap", quietly=TRUE)) {
      stop("plot=TRUE requires package 'pheatmap'")
    }
    if (!requireNamespace("gridExtra", quietly=TRUE)) {
      stop("plot=TRUE requires package 'gridExtra'")
    }
  }
  nsnps <- formatC(sapply(sum_stat,nrow), width=2, flag=0)
  message(paste0("pre:  ",paste(nsnps,collapse=",")))
  gs <- list()
  if (!is.null(score)) {
    stopifnot(is(score, "character"))
    stopifnot(score %in% colnames(sum_stat[[1]]))
  }
  for (j in seq_along(sum_stat)) {
    if (nrow(sum_stat[[j]]) == 1) {
      # need to add the 'collapsed' column anyway
      if (!is.null(snp_id)) {
        sum_stat[[j]]$collapsed <- sum_stat[[j]][[snp_id]]
      }      
      next
    }
    if (plot) {
      gs[[2*j-1]] <- pheatmap::pheatmap(ld_mat[[j]], breaks=seq(-1,1,length=101),
                                        cluster_rows=FALSE, cluster_cols=FALSE,
                                        show_rownames=FALSE, show_colnames=FALSE,
                                        border_color=NA, 
                                        silent=TRUE)$gtable
    }
    hc <- hclust(as.dist(1-abs(ld_mat[[j]])))
    #plot(hc); abline(h=1-threshold,col="red")
    hclusters <- cutree(hc, h=1-threshold)
    nhclust <- max(hclusters)
    # greedy reduction:
    # just take the first SNP for each hierarchical cluster,
    # unless score is specified, then use the SNP with largest score
    if (!is.null(snp_id)) {
      collapsed <- sapply(split(sum_stat[[j]][[snp_id]], hclusters),
                          paste, collapse=",")
    }
    if (!is.null(score)) {
      z <- sum_stat[[j]][[score]]
      tmp <- hclusters
      for (i in seq_len(nhclust)) {
        max.score <- max(z[hclusters == i])
        tmp[hclusters == i][z[hclusters == i] < max.score] <- NA
      }
      hclusters <- tmp
    }
    idx <- match(seq_len(nhclust), hclusters)
    ld_mat[[j]] <- ld_mat[[j]][idx,idx,drop=FALSE]
    sum_stat[[j]] <- sum_stat[[j]][idx,]
    if (!is.null(snp_id)) {
      sum_stat[[j]]$collapsed <- collapsed
    }
    if (two.ld) {
      ld_mat2[[j]] <- ld_mat2[[j]][idx,idx,drop=FALSE]
    }
    if (plot) {
      gs[[2*j]] <- pheatmap::pheatmap(ld_mat[[j]], breaks=seq(-1,1,length=101),
                                      cluster_rows=FALSE, cluster_cols=FALSE,
                                      show_rownames=FALSE, show_colnames=FALSE,
                                      border_color=NA, 
                                      silent=TRUE)$gtable
    }
  }
  if (plot) {
    gridExtra::grid.arrange(ncol=2, grobs=gs)
  }
  nsnps <- formatC(sapply(sum_stat,nrow), width=2, flag=0)
  message(paste0("post: ",paste(nsnps,collapse=",")))
  if (two.ld) {
    list(sum_stat=sum_stat, ld_mat=ld_mat, ld_mat2=ld_mat2)
  } else {
    list(sum_stat=sum_stat, ld_mat=ld_mat)
  }
}

#' Flip alleles and gather results into lists
#'
#' A helper function to flip alleles from eQTL and GWAS
#' datasets, such that they agree on the effect allele,
#' that the SNPs in a signal cluster are in postive
#' correlation with the index SNP (eQTL), and that the
#' effect allele is coded such that it is the expression
#' increasing allele.
#' This is recommended to run after
#' \code{\link{collapseHighCorSNPs}}, and before
#' \code{\link{fitBetaColoc}}.
#'
#' @param sum_stat list of summary statistic tables.
#' A list over signal clusters, where each element is a
#' data.frame with summary statistics from eQTL and GWAS
#' datasets. The names of the columns are specified by
#' arguments below (e.g. \code{a}, \code{b}, \code{ref},
#' \code{eff}, etc.)
#' @param ld_mat list of LD matrices
#' @param ld_mat2 optional second list of LD matrices
#' (for different populations). it will be returned
#' alongside the first \code{ld_mat}, which is used
#' for the allele flipping. the second list of LD
#' matrices is just flipped in the same way
#' @param a name of A in columns of \code{sum_stat} ("eQTL")
#' @param b name of B ("GWAS")
#' @param ref name of reference allele
#' @param eff name of effect allele
#' @param beta name of estimated coefficient
#' @param se name of standard error
#' @param a2_plink name of the column representing the a2 allele
#' (reference allele) according to plink. the default for plink v1.9
#' and earlier was to reset a2 to the major allele, unless
#' an optional flag was used. in plink v2.0 and onward, one
#' should check to see which allele is used as reference for
#' calculating the LD matrix
#' @param a2_plink_mat2 name of the column representing
#' the a2 allele for the second LD matrix, \code{ld_mat2}
#' (needed only if \code{ld_mat2} was specified)
#' @param snp_id name of SNP id
#' @param sep character separator in column names that involve A/B
#' @param ab_last logical, A/B descriptor is last in column names
#' (e.g. "beta_eqtl", "se_eqtl"))
#' @param alleles_same logical, A/B/LD matrix alleles are identical
#' @param plot logical, draw a scatterplot of the flipped betas
#'
#' @return list with estimated coefficients, standard
#' errors, LD matrix, and alleles data.frame
#' 
#' @export
flipAllelesAndGather <- function(sum_stat, ld_mat,
                                 ld_mat2=NULL,
                                 a, b, ref, eff,
                                 beta, se,
                                 a2_plink,
                                 a2_plink_mat2=NULL,
                                 snp_id,
                                 sep, ab_last=TRUE,
                                 alleles_same=FALSE,
                                 plot=TRUE) {

  if (any(duplicated(do.call(rbind, sum_stat)[[snp_id]]))) {
    warning("duplicate SNPs across signal clusters")
  }
  stopifnot(all(sapply(ld_mat, is, "matrix")))
  stopifnot(nrow(ld_mat[[1]]) == ncol(ld_mat[[1]]))
  stopifnot(nrow(ld_mat[[1]]) == nrow(sum_stat[[1]]))
  two.ld <- !is.null(ld_mat2)
  if (two.ld) {
    stopifnot(all(sapply(ld_mat2, is, "matrix")))
    stopifnot(nrow(ld_mat2[[1]]) == ncol(ld_mat2[[1]]))
    stopifnot(nrow(ld_mat2[[1]]) == nrow(sum_stat[[1]]))
    stopifnot(!is.null(a2_plink_mat2))
  }
  # the following allow for arbitrary incoming column names.
  # the point of this is to reduce mistakes that might occur
  # if users manually had to modify their column names.
  if (!alleles_same) {
    ref_a_nm <- if (ab_last) paste(ref, a, sep=sep) else paste(a, ref, sep=sep)
    ref_b_nm <- if (ab_last) paste(ref, b, sep=sep) else paste(b, ref, sep=sep)
    eff_a_nm <- if (ab_last) paste(eff, a, sep=sep) else paste(a, eff, sep=sep)
    eff_b_nm <- if (ab_last) paste(eff, b, sep=sep) else paste(b, eff, sep=sep)
  } else {
    ref_a_nm <- ref
    eff_a_nm <- eff
  }
  beta_a_nm <- if (ab_last) paste(beta, a, sep=sep) else paste(a, beta, sep=sep)
  beta_b_nm <- if (ab_last) paste(beta, b, sep=sep) else paste(b, beta, sep=sep)
  se_b_nm <- if (ab_last) paste(se, b, sep=sep) else paste(b, se, sep=sep)
  se_a_nm <- if (ab_last) paste(se, a, sep=sep) else paste(a, se, sep=sep)
  
  beta_hat_a <- list()
  beta_hat_b <- list()
  se_a <- list()
  se_b <- list()
  Sigma <- list()
  Sigma2 <- list()
  alleles <- list()
  
  for (j in seq_along(sum_stat)) {
    # the index SNP in the signal cluster for A (eQTL) based on Z stat
    idx <- which.max(abs(sum_stat[[j]][[beta_a_nm]]/sum_stat[[j]][[se_a_nm]]))
    if (!alleles_same) {
      # check: reference B (GWAS) allele must be either reference or effect allele in A (eQTL)
      stopifnot(sum_stat[[j]][[ref_b_nm]] == sum_stat[[j]][[ref_a_nm]] |
                sum_stat[[j]][[ref_b_nm]] == sum_stat[[j]][[eff_a_nm]])
      # flip B (GWAS) so that same ref allele is described for B (GWAS) as for A (eQTL)
      flip <- which(sum_stat[[j]][[ref_b_nm]] != sum_stat[[j]][[ref_a_nm]])
      sum_stat[[j]]$beta_b_flipped <- sum_stat[[j]][[beta_b_nm]]
      sum_stat[[j]]$beta_b_flipped[flip] <- -1 * sum_stat[[j]][[beta_b_nm]][flip]
    } else {
      # no flipping, bc alleles the same
      sum_stat[[j]]$beta_b_flipped <- sum_stat[[j]][[beta_b_nm]]
    }
    # flip alleles other than the index so they have positive correlation (LD) with index
    # and based on what the major allele is according to plink
    ld.sign <- sign(ld_mat[[j]][,idx])
    if (!alleles_same) {
      plink.agree <- ifelse(sum_stat[[j]][[a2_plink]] == sum_stat[[j]][[ref_a_nm]], 1, -1)
    } else {
      plink.agree <- rep(1, nrow(sum_stat[[j]]))
    }
    beta_a <- plink.agree * ld.sign * sum_stat[[j]][[beta_a_nm]]
    beta_b <- plink.agree * ld.sign * sum_stat[[j]]$beta_b_flipped
    # only flip LD matrix based on positive correlation with index (bc it comes from plink)
    ld.flipped <- t(t(ld_mat[[j]]) * ld.sign) * ld.sign
    Sigma[[j]] <- ld.flipped
    if (two.ld) {
      two.ld.agree <- ifelse(sum_stat[[j]][[a2_plink]] == sum_stat[[j]][[a2_plink_mat2]], 1, -1)
      ld.sign <- ld.sign * two.ld.agree
      Sigma2[[j]] <- t(t(ld_mat2[[j]]) * ld.sign) * ld.sign
    }
    # record the alleles after all the flipping
    alleles[[j]] <- data.frame(
      id=sum_stat[[j]][[snp_id]],
      ref=sum_stat[[j]][[ref_a_nm]],
      eff=sum_stat[[j]][[eff_a_nm]],
      stringsAsFactors=FALSE)
    # carry over collapsed SNPs if column is present
    if ("collapsed" %in% colnames(sum_stat[[j]])) {
      alleles[[j]]$collapsed <- sum_stat[[j]]$collapsed
    }
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
              Sigma=Sigma)
  if (two.ld) {
    out$Sigma2 <- Sigma2
  }
  out$alleles <- alleles
  if (plot) {
    plotInitEstimates(out, a=a, b=b)
  }
  return(out)
}

#' Plot initial estimates over signal clusters
#'
#' @param x list of signal clusters data with \code{beta_hat_a}
#' and \code{beta_hat_b} lists
#' @param label what preceeds \code{a} and \code{b} in
#' the x- and y-axis labels
#' @param a name of A experiment
#' @param b name of B experiment
#'
#' @export
plotInitEstimates <- function(x, label="Effect size of", a="eQTL", b="GWAS") {
  nsnp <- lengths(x$beta_hat_a)
  plot(unlist(x$beta_hat_a), unlist(x$beta_hat_b),
       xlab=paste(label, a),
       ylab=paste(label, b), 
       col=rep(seq_along(nsnp),nsnp),
       pch=rep(seq_along(nsnp),nsnp))
  text(unlist(x$beta_hat_a), unlist(x$beta_hat_b),
       do.call(c, lapply(nsnp, seq_len)), pos=4, cex=.5,
       col=rep(seq_along(nsnp),nsnp))
  abline(h=0, col=rgb(0,0,0,.3))
}

#' Extract SNPs from colocalization for slope fitting
#'
#' Extracts one or more SNPs from each signal cluster
#' based on the posterior estimate of the effect size
#' for A (largest effect size in the positive direction).
#' 
#' @param res list with the following named elements:
#' \itemize{
#' \item \code{beta_hat_a} - list of point estimates of coefficients for A from colocalization
#' \item \code{beta_hat_b} - " " for B
#' \item \code{sd_a} - list of sampling SD for \code{beta_hat_a} (in practice original
#' SE are provided here)
#' \item \code{sd_b} - " " for \code{beta_hat_b} " "
#' \item {alleles} (optional) list of data.frame with allele information
#' }
#' @param niter number of iterations of EM to run
#' for mclust, if set to 0, only the maximum
#' variant (in terms of A effect size) per
#' signal cluster is output. Default is to not
#' run clustering, but to take the SNP with the
#' largest effect size in A (in the positive direction)
#' @param plot logical, draw a before after of which
#' variants will be included for slope estimation
#' @param label what preceeds \code{a} and \code{b} in
#' the x- and y-axis labels
#' @param a name of A experiment
#' @param b name of B experiment
#'
#' @return list of vectors of the first four arguments,
#' collapsed now across signal clusters, representing
#' variants with positive effect on A. So the null variants
#' have been removed (and any variants per cluster that
#' indicated a negative effect on A). If \code{alleles}
#' data.frames were included in the input, they will
#' also be passed through as a single data.frame with the
#' selected SNPs per signal cluster
#'
#' @export
extractForSlope <- function(res,
                            niter=0,
                            plot=TRUE,
                            label="Effect size of",
                            a="eQTL", b="GWAS") {
  stopifnot(all(c("beta_hat_a","beta_hat_b","sd_a","sd_b") %in% names(res)))
  nsnp <- lengths(res$beta_hat_a)
  stopifnot(all(lengths(res$beta_hat_b) == nsnp))
  stopifnot(all(lengths(res$sd_a) == nsnp))
  stopifnot(all(lengths(res$sd_b) == nsnp))
  if (niter == 0) {
    z <- ifelse(
      unlist(lapply(res$beta_hat_a,
                    function(x) x == max(x))), 2, 1)
  } else {
    if (!requireNamespace("mclust", quietly=TRUE)) {
      stop("niter > 0 requires package 'mclust'")
    }
    beta_max_a <- sapply(res$beta_hat_a, max)
    dat <- pmax(unlist(res$beta_hat_a),0)
    kfit <- kmeans(dat, centers=c(0,mean(beta_max_a)))
    z <- kfit$cluster
    for (i in 1:niter) {
      ms <- mclust::mstepV(data=dat, z=mclust::unmap(z))
      es <- mclust::estepV(data=dat, parameters=ms$parameters)
      z <- ifelse(es$z[,2] > .5, 2, 1)
    }
  }
  if (plot) {
    par(mfrow=c(1,2))
    plot(unlist(res$beta_hat_a), unlist(res$beta_hat_b),
         col=rep(seq_along(nsnp),nsnp),
         pch=rep(seq_along(nsnp),nsnp),
         xlab=paste(label,a), ylab=paste(label,b))
    plot(unlist(res$beta_hat_a), unlist(res$beta_hat_b),
         col=ifelse(z == 1, "black", "blue"),
         pch=rep(seq_along(nsnp),nsnp),
         xlab=paste(label,a), ylab=paste(label,b))
  }
  stopifnot(any(z == 2))
  idx <- z == 2
  out <- list(beta_hat_a=unlist(res$beta_hat_a)[idx],
              beta_hat_b=unlist(res$beta_hat_b)[idx],
              sd_a=unlist(res$sd_a)[idx],
              sd_b=unlist(res$sd_b)[idx])
  if ("alleles" %in% names(res)) {
    alleles <- res$alleles
    alleles <- do.call(rbind, alleles)[idx,,drop=FALSE]
    out$alleles <- alleles
  }
  out
}

#' Make simple simulated summary data
#'
#' @param nsnp number of SNPs per signal cluster
#' @param idx the causal SNP (same per cluster for
#' simplicity)
#' @param alpha the true slope of B coefficients over A coefficients
#' @param sigma the SD of true B coefficients around the conditional
#' values given true A coefficients
#' @param betas the true A coefficients
#' @param se the standard errors for betas
#' @param n_mult how many more samples the B study has
#' 
#' @return a list of \code{beta_hat_a}, \code{beta_hat_b},
#' \code{se_a}, \code{se_b},
#' \code{Sigma_a}, \code{Sigma_b} (themselves lists), and \code{alleles}
#' (a list of data.frames each with
#' \code{id}, \code{ref}, \code{eff} for the SNP id,
#' reference allele, and effect allele).
#'
#' @importFrom MASS mvrnorm
#' @export
makeSimDataForMrlocus <- function(nsnp=c(7:10), idx=5,
                                  alpha=.5, sigma=.05,
                                  betas=1:4, se=.25, n_mult=1) {
  stopifnot(idx >= 3)
  stopifnot(all(nsnp >= idx))
  stopifnot(all(betas > 0))
  stopifnot(length(betas) == length(nsnp))
  nclust <- length(nsnp)
  Sigma_a <- Sigma_b <- list()
  for (j in 1:nclust) {
    Sigma_a[[j]] <- diag(nsnp[j]) # A will be eQTL
    Sigma_b[[j]] <- diag(nsnp[j]) # B will be GWAS
    z <- idx + -2:2
    Sigma_a[[j]][z,z] <- ifelse(Sigma_a[[j]][z,z] == 0, .5, 1)
    Sigma_b[[j]][z,z] <- ifelse(Sigma_b[[j]][z,z] == 0, .5, 1)
  }
  x <- idx - 1
  y <- max(nsnp) - idx
  beta <- lapply(1:nclust, function(j) (rep(c(0,betas[j],0),c(x,1,y)))[1:nsnp[j]])
  beta_hat_a <- beta_hat_b <- beta
  se_a <- se_b <- lapply(1:nclust, function(j) rep(se, nsnp[j]))
  if (n_mult != 1) {
    se_b <- lapply(se_b, function(x) x/sqrt(n_mult))
  }
  mu <- mean(sapply(beta, `[`, x+1))
  for (j in 1:nclust) {
    beta_a_j <- beta[[j]]
    beta_hat_a[[j]] <- MASS::mvrnorm(1,
                       mu=Sigma_a[[j]] %*% beta_a_j,
                       diag(se_a[[j]]) %*% Sigma_a[[j]] %*% diag(se_a[[j]]))
    beta_b_j <- alpha * beta_a_j + ifelse(beta_a_j==0,0,rnorm(nsnp[j],0,sigma))
    beta_hat_b[[j]] <- MASS::mvrnorm(1,
                       mu=Sigma_b[[j]] %*% beta_b_j,
                       diag(se_b[[j]]) %*% Sigma_b[[j]] %*% diag(se_b[[j]]))
  }
  alleles <- lapply(nsnp, function(n) {
    data.frame(id=paste0("rs",round(runif(n,1,1e6))),
               ref=rep(c("A","C"),length.out=n),
               eff=rep(c("T","G"),length.out=n))
  })
  list(beta_hat_a=beta_hat_a,
       beta_hat_b=beta_hat_b,
       se_a=se_a,
       se_b=se_b,
       Sigma_a=Sigma_a,
       Sigma_b=Sigma_b,
       alleles=alleles)
}

#' Plot estimates from MRLocus slope fitting step
#'
#' @param res the output from \code{\link{fitSlope}}
#' @param q the quantiles of the posterior
#' to use for drawing the uncertainty on the slope.
#' The default is an 80 percent interval
#' @param sigma_mult multiplier on estimate of sigma
#' for drawing the dispersion band
#' (e.g. \code{qnorm(1 - .2/2) ~= 1.28} should include
#' 80 percent of coefficient pairs)
#' @param label what preceeds \code{a} and \code{b} in
#' the x- and y-axis labels
#' @param a name of A experiment
#' @param b name of B experiment
#' @param xlim xlim (if NULL will be set automatically)
#' @param ylim ylim (if NULL will be set automatically)
#' @param legend logical, whether to show a legend
#' @param digits number of digits to show in legend
#' @param ... arguments passed to \code{plot}
#'
#' @export
plotMrlocus <- function(res, 
                        q=c(.1,.9),
                        sigma_mult=1.28,
                        label="Effect size of",
                        a="eQTL", b="GWAS",
                        xlim=NULL,
                        ylim=NULL,
                        legend=TRUE,
                        digits=3,
                        ...) {
  stopifnot(length(q) == 2)
  stansum <- rstan::summary(res$stanfit, pars=c("alpha","sigma"), probs=q)$summary
  alpha.hat <- stansum["alpha","mean"]
  qs <- paste0(q * 100,"%")
  alpha.qs <- stansum["alpha",qs]
  sigma.hat <- stansum["sigma","mean"]
  if (is.null(xlim)) {
    xx <- max(res$beta_hat_a)
    xlim <- c(0, 1.5*xx)
  } else {
    xx <- 1.33*xlim[2]
  }
  yy <- 1.5*max(abs(res$beta_hat_b))
  if (is.null(ylim)) {
    ylim <- c(-yy, yy)
  }
  plot(res$beta_hat_a, res$beta_hat_b,
       xlim=xlim, ylim=ylim, type="n",
       xlab=paste(label, a),
       ylab=paste(label, b), ...)

  # alpha (slope), sigma, and uncertainty on slope
  polygon(c(0,2*xx,2*xx,0),
          c(-sigma_mult*sigma.hat, alpha.hat*2*xx -sigma_mult*sigma.hat,
            alpha.hat*2*xx + sigma_mult*sigma.hat, sigma_mult*sigma.hat),
          col=rgb(0,0,1,.1), border=NA)
  segments(0, 0, 2*xx, alpha.hat*2*xx, col="blue", lwd=2)
  segments(0, 0, 2*xx, (alpha.qs[1])*2*xx,
           col=rgb(0,0,1,.5), lwd=2, lty=2)
  segments(0, 0, 2*xx, (alpha.qs[2])*2*xx,
           col=rgb(0,0,1,.5), lwd=2, lty=2)
  abline(h=0, col=rgb(0,0,0,.25))

  # the pairs and their SEs
  points(res$beta_hat_a, res$beta_hat_b, pch=19)
  arrows(res$beta_hat_a - res$sd_a, res$beta_hat_b,
         res$beta_hat_a + res$sd_a, res$beta_hat_b,
         code=3, angle=90, length=.05)
  arrows(res$beta_hat_a, res$beta_hat_b - res$sd_b,
         res$beta_hat_a, res$beta_hat_b + res$sd_b,
         code=3, angle=90, length=.05)

  if (legend) {
    where <- if (alpha.hat > 0) "bottomleft" else "topleft"
    slope.leg <- as.expression(bquote(paste("slope ", hat(alpha)," = ",
                                            .(round(alpha.hat,digits)))))
    slope.int.leg <- paste0(100*diff(q), "% int. = (",
                            round(alpha.qs[1],digits),
                           ", ",round(alpha.qs[2],digits),
                           ")")
    sigma.leg <- if (sigma_mult == 1) {
                   as.expression(bquote(paste(hat(sigma)," = ",
                                              .(round(sigma.hat,digits)))))
                 } else {
                   as.expression(bquote(paste("" %+-% .(sigma_mult) %*% "(", hat(sigma)," = ",
                                              .(round(sigma.hat,digits)), ")")))
                 }
    legend(where,
           lwd=c(2,2,5),
           lty=c(1,2,1),
           col=rgb(0,0,1,c(1,.5,.15)),
           inset=.05,
           y.intersp=1.1,
           bg="white",
           legend=c(slope.leg, slope.int.leg, sigma.leg))
  }
}

#' Basic prior checks on MRLocus slope fit
#'
#' This function provides some basic checks on the
#' strength of the prior in the MRLocus slope fitting
#' Bayesian model. It is not desired that the prior
#' overly influences the posterior inference.
#'
#' The posterior-over-prior SD ratio is calculated
#' and returned in a table, and two plots are made
#' that show parameters drawn from the estimated
#' priors (in MRLocus, priors are estimated from the data).
#' Alternatively, the prior predictive draws themselves can
#' be returned instead of the table (by setting \code{type=2}).
#' 
#' If the posterior-over-prior SD ratio is close to 1
#' for either alpha or sigma, this indicates
#' undesirable influence of the prior on the
#' posterior inference. For comparison,
#' some consider a posterior-prior SD ratio of 0.1 or higher
#' to be described as an 'informative prior'
#' (from Stan wiki on prior choice recommendations). 
#' We note that an 'informative prior' alone is not
#' problematic for MRLocus, and the prior estimation steps
#' have been designed to be informative as
#' to reasonable values for the parameters alpha
#' and sigma.
#'
#' The plots show parameters generated
#' from the prior and the model. The simulated true values of
#' \code{beta_a} and \code{beta_b} are drawn as black
#' circles (summary statistics would then be drawn from
#' these according to the reported SEs, but this step
#' of the model is omitted in this plot).
#' The two plots differ in that the second plot fixed
#' alpha instead of drawing it from the model
#' (so that the prior for sigma can better be visualized).
#' The fitted estimates of \code{beta_a} and \code{beta_b}
#' from the colocalization step are shown as blue X's.
#' One exception where parameters are not drawn from the prior is:
#' \code{beta_a} values are instead drawn as uniform
#' between 0 and 1.1x the maximum value of \code{beta_hat_a}
#' from the colocalization step (for ease of visualization).
#' 
#' @param res output of \code{\link{fitSlope}}
#' @param n integer, for the plot how many data points to simulate
#' @param plot logical, whether to draw the plots
#' @param type integer, return type. by default (\code{type=1})
#' the function returns a table. By setting \code{type=2},
#' the prior predictive draws for alpha, sigma, \code{beta_a},
#' and \code{beta_b} are returned. See Details regarding the
#' simulated draws for \code{beta_a}
#'
#' @return a data.frame with information
#' about prior and posterior SD for alpha and sigma,
#' and two plots are generated (see Details)
#'
#' @export
priorCheck <- function(res, n=200, plot=TRUE, type=1) {
  stopifnot("stanfit" %in% names(res))
  stopifnot("priors" %in% names(res))

  # just to clean code a bit
  z <- res$priors
  # extract fitted slope
  stansum <- rstan::summary(res$stanfit, pars=c("alpha"), probs=c(.1,.9))$summary
  alpha.hat <- stansum["alpha","mean"]
  alpha.ci <- stansum["alpha",c("10%","90%")]
  # sample data from prior and model
  alpha <- rnorm(n, z$mu_alpha, z$sd_alpha)
  sigma <- abs(rnorm(n, 0, z$sd_sigma))
  max_a <- 1.1 * max(res$beta_hat_a)
  beta_a <- runif(n, 0, max_a)
  beta_b <- rnorm(n, alpha * beta_a, sigma)

  # plots:

  if (plot) {
    par(mfrow=c(1,2))
    plot(beta_a, beta_b, main="Prior predictive",
         cex=.75, col=rgb(0,0,0,.5))
    points(res$beta_hat_a, res$beta_hat_b, col="blue", pch=4, lwd=2, cex=2)
    segments(0, 0, max_a, alpha*max_a, col=rgb(0,0,0,.1))
    abline(0, alpha.hat, col="blue", lwd=2)
    abline(0, alpha.ci[1], col=rgb(0,0,1,.5), lwd=2, lty=2)
    abline(0, alpha.ci[2], col=rgb(0,0,1,.5), lwd=2, lty=2)
    abline(h=0, lwd=2)
    
    beta_b_fix_alpha <- rnorm(n, z$mu_alpha * beta_a, sigma)
    plot(beta_a, beta_b_fix_alpha, main="Prior predictive (fixed alpha)",
         cex=.75, col=rgb(0,0,0,.5), ylab="beta_b")
    points(res$beta_hat_a, res$beta_hat_b, col="blue", pch=4, lwd=2, cex=2)
    abline(0, alpha.hat, col="blue", lwd=2)
    abline(0, alpha.ci[1], col=rgb(0,0,1,.5), lwd=2, lty=2)
    abline(0, alpha.ci[2], col=rgb(0,0,1,.5), lwd=2, lty=2)
    abline(h=0, lwd=2)
  }
  
  # table / data.frame with prior draws

  if (type == 1) {
    alpha_post_sd <- rstan::summary(res$stanfit, pars="alpha")$summary[,"sd"]
    sigma_post_sd <- rstan::summary(res$stanfit, pars="sigma")$summary[,"sd"]
    out <- data.frame(
      parameter=c("alpha","sigma"),
      prior_sd=c(res$priors$sd_alpha, res$priors$sd_sigma),
      post_sd=c(alpha_post_sd, sigma_post_sd))
    out$po_pr_ratio <- out[,"post_sd"]/out[,"prior_sd"]
  } else {
    out <- data.frame(alpha=alpha, sigma=sigma, beta_a=beta_a, beta_b=beta_b)
  }
  out
  
}

## addPointers <- function(res, ylim, alpha.hat) {
##   # only works for a pos slope example
##   xr <- 1.5 * max(res$beta_hat_a)
##   yr <- diff(ylim)
##   idx <- which.min(res$beta_hat_a)
##   arrows(res$beta_hat_a[idx] + xr/20, -yr/3,
##          res$beta_hat_a[idx], res$beta_hat_b[idx] - yr/20,
##          angle=45, length=.05)
##   text(res$beta_hat_a[idx] + xr/20, -yr/3,
##        "MRLocus est. coef.\nand SE bars",
##        pos=1, cex=.75)
##   blue <- "blue3"
##   arrows(.5 * xr, -yr/8, .4 * xr,
##          .75 * alpha.hat * .4 * xr,
##          angle=45, length=.05, col=blue)
##   text(.5 * xr, -yr/6,
##        c("80% dispersion\n",
##          as.expression(bquote(paste("band (using ",hat(sigma),")")))),
##        pos=1, cex=.75, col=blue)
##   arrows(.75 * xr, yr/5, .75 * xr,
##          .97 * alpha.hat * .75 * xr,
##          angle=45, length=.05, col=blue)
##   text(.75 * xr, yr/6,
##        c("gene-to-trait\n",
##          as.expression(bquote(paste("slope (",hat(alpha),")")))),
##        pos=1, cex=.75, col=blue)
##   arrows(.9 * xr, -yr/8, .9 * xr,
##          .85 * alpha.hat * .9 * xr,
##          angle=45, length=.05, col=blue)
##   text(.9 * xr, -yr/8,
##        "80% interval\non slope",
##        pos=1, cex=.75, col=blue)
## }
