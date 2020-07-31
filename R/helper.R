#' Collapse high correlation SNPs
#'
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#' @param ld_mat2 optional second list of LD matrices
#' (for different populations). it will be returned
#' alongside the first \code{ld_mat}, which is used
#' for the collapsing. the second list of LD
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
#'
#' @return list with modified ld_mat and sum_stat lists
#' 
#' @export
collapseHighCorSNPs <- function(sum_stat, ld_mat, ld_mat2=NULL,
                                thresh=.95, score=NULL, plot=TRUE) {
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
  nsnps <- formatC(sapply(sum_stat,nrow), width=2, flag=0)
  message(paste0("pre:  ",paste(nsnps,collapse=",")))
  gs <- list()
  if (!is.null(score)) {
    stopifnot(is(score, "character"))
    stopifnot(score %in% colnames(sum_stat[[1]]))
  }
  for (j in seq_along(sum_stat)) {
    if (nrow(sum_stat[[j]]) == 1) next
    if (plot) {
      gs[[2*j-1]] <- pheatmap::pheatmap(ld_mat[[j]], breaks=seq(-1,1,length=101),
                                        cluster_rows=FALSE, cluster_cols=FALSE,
                                        show_rownames=FALSE, show_colnames=FALSE,
                                        border_color=NA, 
                                        silent=TRUE)$gtable
    }
    hc <- hclust(as.dist(1-abs(ld_mat[[j]])))
    #plot(hc); abline(h=1-thresh,col="red")
    hclusters <- cutree(hc, h=1-thresh)
    nhclust <- max(hclusters)
    # greedy reduction:
    # just take the first SNP for each hierarchical cluster,
    # unless score is specified, then use the SNP with largest score
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
#' @param sum_stat list of summary statistic tables
#' @param ld_mat list of LD matrices
#' @param ld_mat2 optional second list of LD matrices
#' (for different populations). it will be returned
#' alongside the first \code{ld_mat}, which is used
#' for the allele flipping. the second list of LD
#' matrices is just flipped in the same way
#' @param a name of A in columns of sum_stat ("eQTL")
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
#' @param snp_id name of SNP id
#' @param sep character separator in column names that involve A/B
#' @param ab_last logical, A/B descriptor is last in column names
#' (e.g. "beta_eqtl", "se_eqtl"))
#' @param alleles_same logical, A/B/LD matrix alleles are identical
#' @param plot logical, draw a scatterplot of the flipped betas
#'
#' @return list with estimated coefficients, standard
#' errors, LD matrix, and allele table
#' 
#' @export
flipAllelesAndGather <- function(sum_stat, ld_mat,
                                 ld_mat2=NULL,
                                 a, b, ref, eff,
                                 beta, se, a2_plink, snp_id,
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
      Sigma2[[j]] <- t(t(ld_mat2[[j]]) * ld.sign) * ld.sign
    }
    # record the alleles after all the flipping
    alleles[[j]] <- data.frame(
      id=sum_stat[[j]][[snp_id]],
      ref=sum_stat[[j]][[ref_a_nm]],
      eff=sum_stat[[j]][[eff_a_nm]],
      stringsAsFactors=FALSE)
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
  if (two.ld) {
    out$Sigma2 <- Sigma2
  }
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
       do.call(c, lapply(nsnp, seq_len)), pos=4, cex=.5,
       col=rep(seq_along(nsnp),nsnp))
  abline(h=0, col=rgb(0,0,0,.3))
}

#' Extract SNPs from colocalization for slope fitting
#'
#' @param res list with the following named elements:
#' \itemize{
#' \item beta_hat_a - vector of length sum(nsnp), first step point estimates of beta for A
#' \item beta_hat_b - " " for B
#' \item sd_a - vector of length sum(nsnp), first step posterior SD (or SE) for beta for A 
#' \item sd_b - " " for B
#' }
#' @param niter number of iterations of EM to run
#' for Mclust, if set to 0, only the maximum
#' variant (in terms of A effect size) per
#' signal cluster is output.
#' @param plot logical, draw a before after of which
#' variants will be included for slope estimation
#' @param a name of A experiment
#' @param b name of B experiment
#'
#' @return list of vectors of the first four arguments,
#' collapsed now across signal clusters, representing
#' variants with positive effect on A. So the null variants
#' have been removed (and any variants per cluster that
#' indicated a negative effect on A)
#'
#' @export
extractForSlope <- function(res,
                            niter=0,
                            plot=TRUE,
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
    beta_max_a <- sapply(res$beta_hat_a, max)
    dat <- pmax(unlist(res$beta_hat_a),0)
    kfit <- kmeans(dat, centers=c(0,mean(beta_max_a)))
    z <- kfit$cluster
    for (i in 1:niter) {
      ms <- mclust::mstep(modelName="V", data=dat, z=unmap(z))
      es <- mclust::estep(modelName="V", data=dat, parameters=ms$parameters)
      z <- ifelse(es$z[,2] > .5, 2, 1)
    }
  }
  if (plot) {
    par(mfrow=c(1,2))
    plot(unlist(res$beta_hat_a), unlist(res$beta_hat_b),
         col=rep(seq_along(nsnp),nsnp),
         pch=rep(seq_along(nsnp),nsnp),
         xlab=paste("beta",a), ylab=paste("beta",b))
    plot(unlist(res$beta_hat_a), unlist(res$beta_hat_b), col=z,
         pch=rep(seq_along(nsnp),nsnp),
         xlab=paste("beta",a), ylab=paste("beta",b))
  }
  stopifnot(any(z == 2))
  idx <- z == 2
  list(beta_hat_a=unlist(res$beta_hat_a)[idx],
       beta_hat_b=unlist(res$beta_hat_b)[idx],
       sd_a=unlist(res$sd_a)[idx],
       sd_b=unlist(res$sd_b)[idx])
}

#' Make simulated data for mrlocus
#'
#' @param nsnp number of SNPs per signal cluster
#' @param idx the causal SNP (same per cluster for
#' simplicity)
#' @param alpha the true slope of B coefficients over A coefficients
#' @param sigma the SD of true B coefficients around the conditional
#' values given true A coefficients
#' @param betas the true A coefficients
#' @param se the standard errors for betas
#' 
#' @return a list of beta_hat_a, beta_hat_b, se_a, and se_b,
#' Sigma_a, and Sigma_b (themselves lists)
#'
#' @export
makeSimDataForMrlocus <- function(nsnp=c(7:10), idx=5,
                                  alpha=.5, sigma=.05,
                                  betas=1:4, se=.25) {
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
  list(beta_hat_a=beta_hat_a,
       beta_hat_b=beta_hat_b,
       se_a=se_a,
       se_b=se_b,
       Sigma_a=Sigma_a,
       Sigma_b=Sigma_b)
}

#' Plot mrlocus estimates
#'
#' @param res the output from fitSlope
#' @param q the quantiles of the posterior
#' to use for drawing the uncertainty on the slope
#' @param a name of A experiment
#' @param b name of B experiment
#' @param ... arguments passed to plot()
#'
#' @export
plotMrlocus <- function(res, 
                        q=c(.1,.9),
                        a="eQTL", b="GWAS",
                        ...) {
  par(mfrow=c(1,1))
  stopifnot(length(q) == 2)
  stansum <- rstan::summary(res$stanfit, pars=c("alpha","sigma"), probs=q)$summary
  alpha.hat <- stansum["alpha","mean"]
  qs <- paste0(q * 100,"%")
  message(paste0("plotting a ",qs[1],"-",qs[2]," interval"))
  alpha.qs <- stansum["alpha",qs]
  sigma.hat <- stansum["sigma","mean"]
  xx <- max(res$beta_hat_a)
  xlim <- c(0, 1.5*xx)
  yy <- 1.5*max(abs(res$beta_hat_b))
  ylim <- c(-yy, yy)
  plot(res$beta_hat_a, res$beta_hat_b,
       xlim=xlim, ylim=ylim, type="n",
       xlab=paste("beta", a),
       ylab=paste("beta", b), ...)

  # alpha (slope), sigma, and uncertainty on slope
  polygon(c(0,2*xx,2*xx,0),
          c(-sigma.hat,alpha.hat*2*xx-sigma.hat,
            alpha.hat*2*xx+sigma.hat,sigma.hat),
          col=rgb(0,0,1,.1), border=NA)
  segments(0, 0, 2*xx, alpha.hat*2*xx, col="blue", lwd=2)
  segments(0, 0, 2*xx, (alpha.qs[1])*2*xx,
           col=rgb(0,0,1,.5))
  segments(0, 0, 2*xx, (alpha.qs[2])*2*xx,
           col=rgb(0,0,1,.5))
  abline(h=0, lty=2)

  # the pairs and their SEs
  points(res$beta_hat_a, res$beta_hat_b, pch=19)
  arrows(res$beta_hat_a - res$sd_a, res$beta_hat_b,
         res$beta_hat_a + res$sd_a, res$beta_hat_b,
         code=3, angle=90, length=.05)
  arrows(res$beta_hat_a, res$beta_hat_b - res$sd_b,
         res$beta_hat_a, res$beta_hat_b + res$sd_b,
         code=3, angle=90, length=.05)
  
  where <- if (alpha.hat > 0) "bottomleft" else "topleft"
  legend(where, lwd=c(2,5), col=rgb(0,0,1,c(1,.15)),
         inset=.05, y.intersp=1.1,
         legend=c(as.expression(bquote(paste(hat(alpha)," = ",
                                             .(round(alpha.hat,3))," [",
                                             .(round(alpha.qs[1],3)),",",
                                             .(round(alpha.qs[2],3)),
                                             "]"))),
                  as.expression(bquote(paste(hat(sigma)," = ",
                                             .(round(sigma.hat,3)))))))
}
