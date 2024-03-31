test_that("penalty works", {
  G <- matrix(c(0, 1, 1, 0, 0,
                1, 0, 1, 0, 0,
                1, 1, 0, 1, 1,
                0, 0, 1, 0, 0,
                0, 0, 1, 0, 1), ncol = 5)

  expect_true(is.numeric(penalty(2, 2, G)))
})

test_that("pseudo-likelihood works", {
  sample_df <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, 5*100, replace = TRUE), ncol = 5)))
  G <- matrix(c(0, 1, 1, 0, 0,
                1, 0, 1, 0, 0,
                1, 1, 0, 1, 1,
                0, 0, 1, 0, 0,
                0, 0, 1, 0, 1), ncol = 5)

  G2 <- matrix(c(0, 1, 1, 1, 0,
                 1, 0, 1, 0, 0,
                 1, 1, 0, 1, 1,
                 1, 0, 1, 0, 0,
                 0, 0, 1, 0, 1), ncol = 5)

  terms <- get_sum_terms(sample_df, G)
  expect_true(is.numeric(sum_terms(terms) + penalty(2, 2, G)))

  updated <- update_sum_terms(terms, sample_df, G2, c(1, 4))

  expect_true(is.numeric(sum_terms(updated) + penalty(2, 2, G2)))
})
