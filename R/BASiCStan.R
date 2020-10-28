#' Stan implementation of BASiCS (with spike-in genes only)
#' @param Data SingleCellExperiment object
#' @param Method Inference method: Variational Bayes ("vb") or Hamiltonian
#' Monte Carlo ("sampling").
#' @param Regression Logical scalar indicating whether or perform the 
#' mean-overdispersion Regression described in Eling et al.
#' @param WithSpikes Do the data contain spike-in genes? See BASiCS for details.
#' @param ... Passed to vb or sampling.
#' @importFrom rstan vb sampling
#' 
#' @export
BASiCStan <- function(
    Data, 
    Method = c("vb", "sampling", "optimizing"), 
    Regression = TRUE,
    WithSpikes = length(altExpNames(Data)) > 0,
    BatchInfo = Data$BatchInfo,
    L = 12,
    ReturnFit = TRUE,
    ...
  ) {

  Method <- match.arg(Method)
  fun <- match.fun(Method)

  if (!is.logical(Regression) && length(Regression) == 1) {
    stop("Invalid value for Regression argument.")
  }
  
  mod <- "basics_regression"

  # if (Regression) {
  #   mod <- paste(mod, "regression", sep = "_")
  # }
  # if (!WithSpikes) {
  #   mod <- paste(mod, "nospikes", sep = "_")
  # }
  model <- stanmodels[[mod]]
  if (!WithSpikes) {
    stop("Not implemented yet")
  }

  if (WithSpikes) {
    spikes <- assay(altExp(Data))
    counts <- counts(Data)
  } else {
    spikes <- NULL
    counts <- counts(Data)
  }
  if (is.null(BatchInfo)) {
    BatchInfo <- rep(1, ncol(Data))
  }

  PP <- BASiCS::BASiCS_PriorParam(
    Data,
    PriorMu = "EmpiricalBayes"
  )
  start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam = PP,
    Regression = TRUE,
    WithSpikes = WithSpikes
  )
  Locations <- BASiCS:::.estimateRBFLocations(start$mu0, L, RBFMinMax = FALSE)

  sdata <- list(
    q = nrow(counts), 
    n = ncol(counts),
    sq = nrow(spikes),
    p = length(unique(BatchInfo)),
    counts = as.matrix(counts),
    spikes = spikes,
    spike_levels = rowData(altExp(Data))[, 2],
    batch_design = model.matrix(~0 + BatchInfo),
    as = 1,
    bs = 1,
    atheta = 1,
    btheta = 1,
    mu_mu = PP$mu.mu,
    smu = sqrt(0.5),
    sdelta = sqrt(0.5),
    aphi = rep(1, ncol(counts)),
    mbeta = rep(0, L),
    ml = Locations[, 1],
    l = L,
    vbeta = diag(L),
    rbf_variance = 1.2,
    eta = 5,
    astwo = 2,
    bstwo = 2
  )
  fit <- fun(model, data = sdata, ...)
  if (ReturnFit) {
    fit
  } else {
    stan2basics(fit, gene_names = rownames(counts), cell_names = colnames(counts))
  }
}

#' Convert stan fits to BASiCS_Chain objects.
#' 
#' @param x A stan object
#' @param gene_names,cell_names Gene and cell names. The reason this argument
#' exists is that by default, stan fit parameters are not named.
#' (NOTE: this should be the same order as the
#' data supplied to BASiCS_stan!!!)
#' 
#' @return A BASiCS_Chain object.
#' @importFrom rstan extract
stan2basics <- function(
    x, 
    gene_names = NULL, 
    cell_names = NULL) {
  
  xe <- extract(x)
  parameters <- list(
    mu = xe[["mu"]],
    delta = xe[["delta"]],
    s = xe[["s"]],
    nu = xe[["nu"]],
    theta = as.matrix(xe[["theta"]])
  )
  for (param in c("epsilon", "phi", "beta")) {
    if (!is.null(xe[[param]])) {
      parameters[[param]] <- xe[[param]]
    }
  }
  gp <- intersect(c("mu", "delta", "epsilon"), names(parameters))
  cp <- intersect(c("s", "nu", "phi"), names(parameters))
  if (!is.null(xe[["beta"]])) {
    parameters[["beta"]] <- xe[["beta"]]
  }
  if (is.null(gene_names)) {
    gene_names <- paste("Gene", seq_len(ncol(xe[["mu"]])))
  }
  if (is.null(cell_names)) {
    cell_names <- paste("Cell", seq_len(ncol(xe[["nu"]])))
  }
  parameters[gp] <- lapply(gp, function(x) {
    colnames(parameters[[x]]) <- gene_names
    parameters[[x]]
  })
  parameters[cp] <- lapply(cp, function(x) {
    colnames(parameters[[x]]) <- cell_names
    parameters[[x]]
  })

  new("BASiCS_Chain", parameters = parameters)
}


vb <- function(..., tol_rel_obj = 1e-3) {
  rstan::vb(..., tol_rel_obj = tol_rel_obj)
}
