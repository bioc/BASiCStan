#' The 'BASiCStan' package.
#'
#' @description Provides an interface to infer the parameters of BASiCS using the
#'  variational inference (ADVI), Markov chain Monte Carlo (NUTS),
#'  and maximum a posteriori (BFGS) inference engines in the Stan programming
#'  language. BASiCS is a Bayesian hierarchical model that uses an adaptive
#'  Metropolis within Gibbs sampling scheme. Alternative inference methods
#'  provided by Stan may be preferable in some situations, for example for
#'  particularly large data or posterior distributions with difficult
#'  geometries. See also \link[BASiCS]{BASiCS_MCMC}
#'
#' @docType package
#' @name BASiCStan-package
#' @useDynLib BASiCStan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom SingleCellExperiment counts logcounts sizeFactors altExp
#' altExpNames
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stats model.matrix setNames
#' @importFrom BASiCS BASiCS_PriorParam
#' @importFrom utils capture.output
#' @importFrom scran calculateSumFactors
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstantools posterior_predict
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan.
#' R package version 2.21.2. \url{https://mc-stan.org}
#'
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#' \url{https://doi.org/10.1371/journal.pcbi.1004333}.
#' 
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
#' \url{https://doi.org/10.1186/s13059-016-0930-3}
#' 
#' Eling et al (2018). Cell Systems.
#' \url{https://doi.org/10.1016/j.cels.2018.06.011}
NULL
