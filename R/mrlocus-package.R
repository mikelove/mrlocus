#' MRLocus - see this page for typical order of functions
#'
#' @description Mendelian Randomization per locus, leveraging eQTL
#' and GWAS summary statistics, for estimation of gene-to-trait
#' effect size and dispersion.
#'
#' The main functions (in order of typical usage) are:
#'
#' \itemize{
#' \item \code{\link{collapseHighCorSNPs}} - collapse high correlation SNPs
#' \item \code{\link{flipAllelesAndGather}} - flip alleles and gather for colocalization
#' \item \code{\link{fitBetaColoc}} - perform colocalization across signal clusters
#' \item \code{\link{extractForSlope}} - extract one SNP per signal cluster for MR analysis
#' \item \code{\link{fitSlope}} - perform MR analysis to estimate gene-to-trait effect
#' \item \code{\link{plotMrlocus}} - plot estimates from MR analysis
#' \item \code{\link{priorCheck}} - perform prior predictive checks
#' }
#'
#' @docType package
#' @name mrlocus-package
#' @aliases mrlocus
#' @useDynLib mrlocus, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#'
#' Anqi Zhu*, Nana Matoba*, Emmaleigh Wilson, Amanda L. Tapia, Yun Li,
#' Joseph G. Ibrahim, Jason L. Stein, Michael I. Love.
#' MRLocus: identifying causal genes mediating a trait through
#' Bayesian estimation of allelic heterogeneity. (2020) bioRxiv
#' \url{https://doi.org/10.1101/2020.08.14.250720}
#' 
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL
