
#' get_symmetric_matrix_from_triangle
#'
#' Constructs a symmetric matrix from a triangle vector.
#'
#' @param tri A vector representing the lower triangle of the matrix.
#' @param d The number of vertices.
#'
#' @return Returns a symmetric matrix constructed from the triangle vector.
#'
#' @export
#'
#' @examples
#' # Example usage
#' triangle <- c(1, 2, 3, 4, 5, 6)  # Example triangle vector
#' d <- 4                            # Number of vertices
#'
#' # Get the symmetric matrix from the triangle
#' symmetric_matrix <- get_symmetric_matrix_from_triangle(triangle, d)
#' print(symmetric_matrix)
get_symmetric_matrix_from_triangle <- function(tri, d){
  tri <- as.numeric(tri)
  m <- matrix(0, ncol = d, nrow = d)
  m[lower.tri(m)] <- tri
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  m
}

#' Exact algorithm that performs the exhaustive algorithm
#'
#' @param sample_df sample data
#' @param card_A cardinality of the alphabet
#' @param d number of vertices
#' @param lambda penalty term
#'
#' @return
#' @export
#' @import tidyverse
#'
exact <- function(sample_df, card_A, d, lambda){
  exact_parallel(sample_df, card_A, d, lambda)
}


#' Exact algorithm that performs the exhaustive algorithm in parallel
#'
#' This function implements an exact algorithm for graph estimation using parallel computing.
#'
#' @param sample_df Sample data.
#' @param card_A Cardinality of the alphabet.
#' @param d Number of vertices.
#' @param lambda Penalty term.
#'
#' @return A list containing the estimated adjacency matrix (G_hat) and the corresponding penalized pseudolikelihood value (logLG_hat).
#'
#' @export
#'
#' @examples
#' # Example usage
#' sample_data <- data.frame(V1 = c(0, 1, 0, 1), V2 = c(1, 0, 1, 0), V3 = c(1, 0, 0, 1))
#' card_A <- 2
#' d <- 3
#' lambda <- 0.5
# exact_parallel(sample_data, card_A, d, lambda)
#'
exact_parallel <- function(sample_df, card_A, d, lambda){
  future::plan(future::multisession, workers = 3)

  vec <- expand.grid(rep(list(0:(card_A-1)),
                         (d*d-d)/2)) # number of terms in the upper triangle

  vec <- tibble::tibble(vec = asplit(vec, 1)) %>%
    dplyr::mutate(tri = purrr::pmap(list(x = vec), ~ as.numeric(..1)))

  result_parallel <- foreach::foreach(i = 1:nrow(vec), .inorder = FALSE,
                                      .options.future = list(packages=c("foreach", "doParallel",
                                                                        "magrittr", "MixingGraph"))
  ) %dofuture% {

    G <- get_synmetric_matrix_from_triangle(vec$tri[[i]], d)
    logLG <- get_and_sum_terms_with_penalty(sample_data = sample_df,
                                            G = G,
                                            card_A = card_A,
                                            lambda = lambda)
    tibble::tibble(logLG = logLG, G = list(G))
  }

  res <- tibble::enframe(result_parallel) %>%
    dplyr::select(-name) %>%
    tidyr::unnest(value)

  G_hat <- (res %>%
               dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
               dplyr::pull(G))[[1]]

  logLG_hat <- (res %>%
                   dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
                   dplyr::pull(logLG))[1]

  future::resetWorkers(future::plan(future::sequential))

  tibble::tibble(G_hat = list(G_hat), logLG_hat = logLG_hat)
}

iterate_possible_graphs <- function(d, sample_df) {
  max_graphs <- 2^(d*(d-1)/2)

  for (i in 1:max_graphs) {
    cat(i, "of", max_graphs, "\n")
    graph <- matrix(0, nrow = d, ncol = d)
    binary_rep <- as.integer(intToBits(i))[1:d*(d-1)/2]
    graph[lower.tri(graph)] <- binary_rep
    graph <- graph + t(graph)  # Make the matrix symmetric

    if(i == 1){
      logLG_hat <- get_and_sum_terms(data$data[3][[1]], graph) - penalty(1.5, 2, graph)
    }else{
      logLG_tmp <- get_and_sum_terms(data$data[3][[1]], graph) - penalty(1.5, 2, graph)

      if(logLG_tmp > logLG_hat){
        G_hat <- graph
        logLG_hat <- logLG_tmp
      }
    }
  }

  graphs
}

iterate_possible_graphs_parallel <- function(n, num_cores = NULL) {
  if (is.null(num_cores)) {
    num_cores <- detectCores()-2
  }

  registerDoParallel(num_cores)

  graphs <- foreach(i = 1:2^(n*(n-1)/2), .combine = "c") %dopar% {
    graph <- matrix(0, nrow = n, ncol = n)
    binary_rep <- as.integer(intToBits(i))[1:n*(n-1)/2]
    graph[lower.tri(graph)] <- binary_rep
    graph <- graph + t(graph)  # Make the matrix symmetric
    graph
  }

  stopImplicitCluster()

  graphs
}
