test_that("greedy: remove and add edge works", {
  G_complete <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
  G_empty <- matrix(0, nrow = 5, ncol = 5)

  G1 <- remove_edge(G_complete, 1,2)

  G2 <- add_edge(G_empty, 1,2)

  expect_true(G1[1,2] == 0 & G1[2,1] == 0)
  expect_true(G2[1,2] == 1 & G2[2,1] == 1)
})

test_that("greedy: backwards works", {
  set.seed(1)
  d <- 3
  sample_df <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, d*20, replace = TRUE),
                                                        ncol = d)))

  res <- greedy_backwards(sample_df, lambda = 2, card_A = 2, d = d)

  expect_true(is.matrix(res$G_hat[[1]]))
  expect_true(is.numeric(res$logLG_hat))
})

test_that("greedy: forward works", {
  set.seed(1)
  d <- 3
  sample_df <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, d*20, replace = TRUE),
                                                        ncol = d)))

  res <- greedy_forward(sample_df, lambda = 1.5, card_A = 2, d = d)

  expect_true(is.matrix(res$G_hat[[1]]))
  expect_true(is.numeric(res$logLG_hat))
})

test_that("greedy: forward all results works", {
  set.seed(1)
  d <- 3
  sample_df <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, d*20, replace = TRUE),
                                                        ncol = d)))

  res <- greedy_forward(sample_df, lambda = 2, card_A = 0, d = d, TRUE)

  expect_true(nrow(res) > 1)
})
