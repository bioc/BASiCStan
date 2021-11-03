set.seed(42)
mock <- BASiCS::BASiCS_MockSCE(NGenes = 5, NCells = 5, NSpikes = 5)

test_that("stan works in basic modes", {
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
