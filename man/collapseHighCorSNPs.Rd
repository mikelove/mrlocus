% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{collapseHighCorSNPs}
\alias{collapseHighCorSNPs}
\title{Collapse high correlation SNPs}
\usage{
collapseHighCorSNPs(
  sum_stat,
  ld_mat,
  ld_mat2 = NULL,
  threshold = 0.95,
  score = NULL,
  plot = TRUE,
  snp_id = NULL
)
}
\arguments{
\item{sum_stat}{list of summary statistic tables,
which is a list over signal clusters.
Each element of the list should be a data.frame
describing the eQTL and GWAS summary statistics.
The only column in \code{sum_stat} that is used
by the function is \code{score} (optional)}

\item{ld_mat}{list of LD matrices across
signal clusters}

\item{ld_mat2}{optional second list of LD matrices
(for different populations). it will be returned
alongside the first \code{ld_mat}, which is used
for the collapsing. The second list of LD
matrices is just subset to the same set of SNPs
as the first}

\item{threshold}{threshold on absolute value of
correlation for collapsing, e.g. will collapse
if SNPs are more correlated (or anti-correlated)
than this amount}

\item{score}{name of a column of sum_stat data.frames
with a score, such that collapsing will choose the
highest score SNP per collapsed cluster. Otherwise,
if set to NULL, the first SNP will be used}

\item{plot}{logical, draw a before/after grid of plots}

\item{snp_id}{name of SNP id column in \code{sum_stat},
if specified will output a column \code{collapsed}
that lists which SNP ids are represented in the output
(i.e. which highly correlated SNPs were collapsed).}
}
\value{
list with subset \code{ld_mat} and \code{sum_stat}
lists (and \code{ld_mat2} if provided)
}
\description{
A helper function to collapse sets of highly correlated
SNPs within signal clusters. This is recommended to run
before \code{\link{flipAllelesAndGather}}, and before
\code{\link{fitBetaColoc}}.
}
