#' The 'BASiCStan' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
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
#' R package version 2.21.2. https://mc-stan.org
#'
NULL
