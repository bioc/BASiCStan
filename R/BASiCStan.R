#' Stan implementation of BASiCS.
#'
#' The stan programming language enables the use of MAP, VB, and HMC inference.
#' Only the regression mode featuring a joint prior between mean and
#' overdispersion parameters is implemented
#'
#' @param Data SingleCellExperiment object
#' @param Method Inference method. One of: \code{"vb"} for Variational Bayes,
#' \code{"sampling"} for Hamiltonian Monte Carlo,
#' \code{"optimizing"} for or maximum \emph{a posteriori} estimation.
#' @param WithSpikes Do the data contain spike-in genes? See BASiCS for details.
#' When \code{WithSpikes=FALSE}, the cell-specific scaling normalisation 
#' factors are fixed; use \code{NormFactorFun} to specify how size factors
#' should be generated or extracted.
#' @param Regression Use joint prior for mean and overdispersion parameters?
#' Included for compatibility with \code{\link[BASiCS]{BASiCS_MCMC}},
#' but only \code{TRUE} is supported.
#' @param BatchInfo Vector describing which batch each cell is from.
#' @param L Number of regression terms (including slope and intercept) to use in
#' joint prior for mu and delta.
#' @param PriorMu Type of prior to use for mean expression. Default is
#' "EmpiricalBayes", but "uninformative" is the prior used in Eling et al. and
#' previous work.
#' @param NormFactorFun Function that returns cell-specific scaling
#' normalisation factors. See \code{\link[scran]{computeSumFactors}} for
#' details on the default.
#' @param ReturnBASiCS Should the object be converted into a
#' \linkS4class{BASiCS_Chain} object?
#' @param Verbose Should output of the stan commands be printed to
#' the terminal?
#' @param ... Passed to vb or sampling.
#' @importFrom rstan vb sampling
#'
#' @return An object of class \code{\linkS4class{BASiCS_Chain}}.
#' @examples
#' library("BASiCS")
#' sce <- BASiCS_MockSCE(NGenes = 10, NCells = 10)
#'
#' fit_spikes <- BASiCStan(sce)
#'
#' ## uses fixed scaling normalisation factors
#' fit_nospikes <- BASiCStan(sce, WithSpikes = FALSE)
#'
#' @importFrom BASiCS BASiCS_PriorParam
#' @importFrom utils capture.output
#' @export
BASiCStan <- function(
    Data,
    Method = c("vb", "sampling", "optimizing"),
    WithSpikes = length(altExpNames(Data)) > 0,
    Regression = TRUE,
    BatchInfo = Data$BatchInfo,
    L = 12,
    PriorMu = c("EmpiricalBayes", "uninformative"),
    NormFactorFun = scran::calculateSumFactors,
    ReturnBASiCS = TRUE,
    Verbose = TRUE,
    ...
  ) {

    if (!Regression) {
        stop("Only Regression=TRUE is supported.")
    }
    Method <- match.arg(Method)
    fun <- get(Method, mode = "function")
    PriorMu <- match.arg(PriorMu)

    mod <- "basics_regression"

    if (!WithSpikes) {
        mod <- paste(mod, "nospikes", sep = "_")
    }
    model <- stanmodels[[mod]]
    # if (!WithSpikes) {
    #     stop("Not implemented yet")
    # }

    if (WithSpikes) {
        spikes <- assay(altExp(Data))
        counts <- counts(Data)
    } else {
        spikes <- NULL
        counts <- counts(Data)
    }
    if (is.null(BatchInfo) | length(unique(BatchInfo)) == 1) {
        BatchInfo <- 1
        batch_design <- matrix(1, nrow = ncol(Data))
    } else {
        batch_design <- model.matrix(~0 + factor(BatchInfo))
    }

    PP <- BASiCS_PriorParam(
        Data,
        PriorMu = "EmpiricalBayes"
    )
    start <- BASiCS:::.BASiCS_MCMC_Start(
        Data,
        PriorParam = PP,
        Regression = TRUE,
        WithSpikes = WithSpikes
    )
    if (PriorMu == "EmpiricalBayes") {
        PP$mu.mu <- BASiCS:::.EmpiricalBayesMu(Data, 0.5, WithSpikes)
    }
    Locations <- BASiCS:::.estimateRBFLocations(start$mu0, L, RBFMinMax = FALSE)
    size_factors <- match.fun(NormFactorFun)(Data)

    sdata <- list(
        q = nrow(counts),
        n = ncol(counts),
        p = length(unique(BatchInfo)),
        counts = as.matrix(counts),
        batch_design = batch_design,
        size_factors = size_factors,
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
    if (WithSpikes) {
        sdata <- c(
            list(
                sq = nrow(spikes),
                spikes = spikes,
                spike_levels = rowData(altExp(Data))[, 2]
            ),
            sdata
        )
    }
    if (Verbose) {
        fit <- fun(
            model,
            data = sdata,
            init = .init(Data, L, length(unique(BatchInfo)), size_factors),
            ...
        )
    } else {
        capture.output(
            fit <- fun(
                model,
                data = sdata,
                init = .init(Data, L, length(unique(BatchInfo)), size_factors),
                ...
            )
        )
    }
    if (ReturnBASiCS) {
        .stan2basics(
            fit,
            gene_names = rownames(counts),
            cell_names = colnames(counts),
            size_factors = size_factors
        )
    } else {
        fit
    }
}

#' Convert stan fits to BASiCS_Chain objects.
#'
#' @param x A stan object
#' @param gene_names,cell_names Gene and cell names. The reason this argument
#' exists is that by default, stan fit parameters are not named.
#' NOTE: this must be the same order as the
#' data supplied to BASiCS_stan.
#' @param size_factors Cell-specific scaling normalisation factors, to be
#' stored as part of the chain object when \code{WithSpikes=FALSE}.
#'
#' @return A BASiCS_Chain object.
#' @importFrom rstan extract
.stan2basics <- function(
    x,
    gene_names = NULL,
    cell_names = NULL,
    size_factors = NULL) {

    xe <- extract(x)
    parameters <- list(
        mu = xe[["mu"]],
        delta = xe[["delta"]],
        s = xe[["s"]],
        nu = xe[["nu"]],
        theta = xe[["theta"]]
    )
    if (is.null(parameters$nu)) {
        parameters$nu <- t(replicate(nrow(parameters$mu), size_factors))
        parameters$theta <- matrix(1, nrow = nrow(parameters$mu), ncol = 1)
        parameters$s <- parameters$nu
    }
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
    parameters[cp] <- lapply(cp,
        function(x) {
            colnames(parameters[[x]]) <- cell_names
            parameters[[x]]
        }
    )

    new("BASiCS_Chain", parameters = parameters)
}

.init <- function(Data, L, P, SizeFactors) {
    Data <- scuttle::logNormCounts(Data, size.factors = SizeFactors)
    lc <- logcounts(Data)
    gp <- glmGamPoi::glm_gp(Data)
    list(
        mu = exp(rowMeans(lc)),
        delta = gp$overdispersions + 1e-2,
        epsilon = rep(0, nrow(Data)),
        theta = as.array(rep(1, P)),
        beta = rep(0, L),
        lambda = rep(1, nrow(Data)),
        stwo = 1,
        nu = SizeFactors,
        s = rep(1, ncol(Data)),
        phi = rep(1, ncol(Data))
    )
}

vb <- function(..., tol_rel_obj = 1e-3) {
    rstan::vb(..., tol_rel_obj = tol_rel_obj)
}

sampling <- function(...) {
    rstan::sampling(...)
}

optimizing <- function(...) {
    rstan::optimizing(..., verbose = TRUE)
}

extract <- function(fit) {
    UseMethod("extract")
}

#' @export
extract.stanfit <- function(fit) {
    rstan::extract(fit)
}

#' @export
extract.list <- function(fit) {
    pars <- gsub("\\[(\\d+,)*\\d+\\]", "", names(fit$par))
    dims <- gsub("^[a-z0_]+\\[(.*)\\]", "\\1", names(fit$par))
    out <- lapply(unique(pars),
        function(p) {
            ind <- pars == p
            if (!all(grepl("[", names(fit$par)[ind], fixed = TRUE))) {
                return(matrix(fit$par[ind], ncol = 1))
            }
            l <- strsplit(dims[ind], ",")
            d <- do.call(rbind, l)
            d <- matrix(as.numeric(d), ncol = ncol(d))
            array(
                fit$par[ind],
                dim = c(1, apply(d, 2, max))
            )
        }
    )
    setNames(out, unique(pars))
}
