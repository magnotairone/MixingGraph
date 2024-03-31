test_that("SA - neighbor: ", {
  G_complete <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
  G_empty <- matrix(0, nrow = 5, ncol = 5)


  set.seed(1)
  expect_true(abs(sum(G_complete)/2 - sum(choose_neighbor(G_complete))/2) <= 1)
  expect_true(abs(sum(G_empty)/2 - sum(choose_neighbor(G_empty))/2) <= 1)
})

test_that("SA - neighbor: ", {
  G_complete <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
  G_empty <- matrix(0, nrow = 5, ncol = 5)

  set.seed(1)
  res <- identify_diff(G_complete, choose_neighbor(G_complete))
  expect_true(length(res) == 2)
})

test_that("SA - algorithm: ", {
  n_nodes <- 3
  G_0 <- matrix(1, nrow = n_nodes, ncol = n_nodes) - diag(1, n_nodes)
  card_A <- 2
  N <- 100
  C <- 1
  lambda <- 2
  sample_data <- suppressWarnings(tidyr::as_tibble(matrix(sample(0:1, n_nodes*20, replace = TRUE), ncol = n_nodes)))

  r <- SimAn(sample_data, lambda, card_A, G_0, N, C,
        include_logLG = FALSE)

  expect_true(sum(dim(r) == c(n_nodes, n_nodes)) == 2)
})
