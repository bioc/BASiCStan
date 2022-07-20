set.seed(31)
mock <- BASiCS::BASiCS_MockSCE(NGenes = 20, NCells = 20, NSpikes = 10)

test_that("smoke test", {
    expect_s4_class(
        suppressWarnings(BASiCStan(mock, tol_rel_obj = 1, Verbose = FALSE)),
        "BASiCS_Chain"
    )
    expect_s4_class(
        suppressWarnings(
            BASiCStan(
                mock,
                WithSpikes = FALSE, tol_rel_obj = 1, Verbose = FALSE
        )
        ),
        "BASiCS_Chain"
    )
    expect_s4_class(
        suppressWarnings(
            BASiCStan(mock, Method = "optimizing", Verbose = FALSE)
        ),
        "BASiCS_Chain"
    )
    expect_s4_class(
        suppressWarnings(
            BASiCStan(
                mock,
                WithSpikes = FALSE, Method = "optimizing", Verbose = FALSE
            )
        ),
        "BASiCS_Chain"
    )
    expect_s4_class(
        suppressWarnings(
            BASiCStan(
                mock,
                Method = "sampling", iter = 200, chains = 1, Verbose = FALSE
            )
        ),
        "BASiCS_Chain"
    )
    expect_s4_class(
        suppressWarnings(
            BASiCStan(
                mock,
                WithSpikes = FALSE, Method = "sampling",
                iter = 200, chains = 1, Verbose = FALSE
            )
        ),
        "BASiCS_Chain"
    )

})

set.seed(66)
mock <- BASiCS::BASiCS_MockSCE(NGenes = 20, NCells = 20, NSpikes = 10)

test_that("accuracy test", {
    set.seed(42)
    a <- suppressMessages(capture.output(
        fitb <- BASiCS_MCMC(
            mock, N = 20000, Thin = 10, Burn = 10000, Regression = TRUE
        )
    ))
    fit1 <- suppressWarnings(
        BASiCStan(mock, tol_rel_obj = 1e-3, Verbose = FALSE)
    )
    fit2 <- suppressWarnings(
        BASiCStan(mock, WithSpikes = FALSE, tol_rel_obj = 1e-3, Verbose = FALSE)
    )
    expect_equal(
        colMedians(fitb@parameters$mu), colMedians(fit1@parameters$mu),
        tolerance = 1e-1
    )
    expect_gt(
        cor(
            colMedians(fitb@parameters$delta),
            colMedians(fit1@parameters$delta),
            use = "complete.obs"
        ),
        0.7
    )
    expect_gt(
        cor(
            colMedians(fitb@parameters$epsilon),
            colMedians(fit1@parameters$epsilon),
            use = "complete.obs"
        ),
        0.6
    )
    ## no spikes needs offset correcting
    fit2 <- BASiCS:::.offset_correct(fit2, fitb)
    expect_equal(
        colMedians(fitb@parameters$mu), colMedians(fit2@parameters$mu),
        tolerance = 1e-1
    )
    expect_gt(
        cor(
            colMedians(fitb@parameters$delta),
            colMedians(fit2@parameters$delta),
            use = "complete.obs"
        ),
        0.7
    )
    expect_gt(
        cor(
            colMedians(fitb@parameters$epsilon),
            colMedians(fit2@parameters$epsilon),
            use = "complete.obs"
        ),
        0.45
    )
})
