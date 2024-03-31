test_that("exact works", {
  set.seed(1)
  d <- 3
  sample_df <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, d*20, replace = TRUE),
                                                        ncol = d)))

  res <- exact(sample_df, lambda = 2, card_A = 2, d = d)

  expect_true(is.matrix(res$G_hat[[1]]))
  expect_true(is.numeric(res$logLG_hat))
})

exact
