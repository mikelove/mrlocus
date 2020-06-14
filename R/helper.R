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
