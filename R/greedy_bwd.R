#' Remove edge from an adjacency matrix
#'
#' This function removes the edge between two vertices in an adjacency matrix.
#'
#' @param G Adjacency matrix.
#' @param v1 Vertex 1.
#' @param v2 Vertex 2.
#'
#' @return An adjacency matrix with no edge between v1 and v2.
#' @export
#'
#' @examples
#' # Example usage
#' G_complete <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
#' remove_edge(G_complete, 1, 2)
#'
remove_edge <- function(G, v1, v2){
  G[v1, v2] <- 0
  G[v2, v1] <- 0
  G
}

#' List existing edges of an adjacency matrix
#'
#' This function lists all existing edges in an adjacency matrix.
#'
#' @param G Adjacency matrix of a graph.
#'
#' @return A tibble with a list of all edges in graph G.
#' @export
#'
#' @examples
#' # Example usage
#' G <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
#' list_existing_edges(G)
#'
list_existing_edges <- function(G){
  suppressWarnings(
    G %>%
      tibble::as_tibble() %>%
      tibble::rownames_to_column("a") %>%
      tidyr::pivot_longer(names_to = "b",
                          values_to = "edge", cols = -a) %>%
      dplyr::mutate(b = as.integer(stringr::str_remove(b, "V")),
                    a = as.integer(a)) %>%
      dplyr::filter(edge == 1, a < b) %>%
      dplyr::select(-edge)
  )
}

#' Step down of the backwards stepwise greedy algorithm to estimate the graph
#'
#' This function performs a step down in the backwards stepwise greedy algorithm to estimate the graph.
# It removes one edge from the adjacency matrix of the graph and evaluates the change in the log-likelihood.
#'
#' @param sample_df Sample data. The columns should be named as V1, V2, ...
#' @param G Adjacency matrix.
#' @param card_A Cardinality of the alphabet.
#' @param lambda Term in the penalty.
#'
#' @return A tibble containing the adjacency matrix with the removed edge and its corresponding log-likelihood.
#' @export
#'
#' @examples
#' # Example usage
#' sample_data <- data.frame(V1 = c(1, 2, 3), V2 = c(2, 3, 1))
#' G <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
#' card_A <- 3
#' lambda <- 0.01
#' step_down(sample_data, G, card_A, lambda)
#'
step_down <- function(sample_df, G, card_A, lambda) {

  res <- list_existing_edges(G) %>%
    dplyr::bind_cols(tibble::tibble(card_A = card_A,
                                    lambda = lambda)) %>%
    dplyr::mutate(G_minus = purrr::pmap(list(G = list(G), v1 = a, v2 = b), remove_edge),
                  logLG = purrr::pmap(list(sample_data = list(sample_df),
                                           G = G_minus,
                                           card_A = card_A,
                                           lambda = lambda),
                                      get_and_sum_terms_with_penalty)) %>%
    tidyr::unnest(logLG) %>%
    dplyr::select(logLG, G_minus)

  res_max <- res %>%
    dplyr::slice_max(logLG, n = 1, with_ties = FALSE)

  G_minus <- (res_max %>%
                dplyr::pull(G_minus))[[1]]

  logLG_minus <- res_max %>%
    dplyr::pull(logLG)

  return(tibble::tibble(G = list(G_minus), logLG = logLG_minus))
}

#' Step down of the backwards stepwise greedy algorithm to estimate the graph in parallel
#'
#' This function performs a step down in the backwards stepwise greedy algorithm to estimate the graph.
# It removes one edge from the adjacency matrix of the graph and evaluates the change in the log-likelihood.
# This version of the function uses parallel computation for faster execution.
#'
#' @param sample_df Sample data.
#' @param G Adjacency matrix.
#' @param card_A Cardinality of the alphabet.
#' @param lambda Term in the penalty.
#' @param clusters Number of clusters for parallel computation.
#'
#' @return A tibble containing the adjacency matrix with the removed edge and its corresponding log-likelihood.
#' @export
step_down_parallel <- function(sample_df, G, card_A, lambda, clusters = 3) {
  future::plan(future::multisession, workers = clusters)

  possible_edges <- (list_existing_edges(G))

  result_parallel <- foreach::foreach(i = 1:nrow(possible_edges), .inorder = FALSE,
                        .options.future = list(packages=c("foreach", "doParallel", "magrittr", "MixingGraph"))
  ) %dofuture% {
    G_minus <- remove_edge(G, v1 = possible_edges$a[i], v2 = possible_edges$b[i])

    logLG <- get_and_sum_terms_with_penalty(sample_data = sample_df,
                                            G = G_minus,
                                            card_A = card_A,
                                            lambda = lambda)

    tibble::tibble(logLG = logLG, G_plus = list(G_plus))
  }

  res <- tibble::enframe(result_parallel) %>%
    dplyr::select(-name) %>%
    tidyr::unnest(value)

  G_minus <- (res %>%
               dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
               dplyr::pull(G_plus))[[1]]


  logLG_minus <- (res %>%
                   dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
                   dplyr::pull(logLG))[1]

  future::resetWorkers(future::plan(future::sequential))

  return(tibble::tibble(G = list(G_minus), logLG = logLG_minus))
}

#' Implements the greedy backwards stepwise selection algorithm
#'
#' This function implements the greedy backwards stepwise selection algorithm to estimate the graph structure.
#' It iteratively removes one edge at a time from the initial graph until no improvement in log-likelihood is observed.
#'
#' @param sample_df Sample data.
#' @param lambda Penalty term.
#' @param card_A Cardinality of the alphabet.
#' @param d Number of vertices in the graph.
#' @param show_all_steps Logical, indicating whether to show all steps of the algorithm.
#'
#' @return A tibble containing the estimated adjacency matrix and its corresponding log-likelihood.
#' @export
#' @examples
#' # Example usage
#' sample_data <- data.frame(V1 = c(1, 2, 3), V2 = c(2, 3, 1))
#' lambda <- 0.01
#' card_A <- 3
#' d <- 2
#' show_all_steps <- TRUE
#' greedy_backwards(sample_data, lambda, card_A, d, show_all_steps)
#'
greedy_backwards <- function(sample_df, lambda, card_A, d, show_all_steps = FALSE) {

  G_hat <- matrix(1, nrow = d, ncol = d) - diag(1, d)
  logLG_hat <- get_and_sum_terms_with_penalty(sample_df, G_hat, card_A, lambda)

  G_complete <- matrix(1, nrow = d, ncol = d) - diag(1, d)

  G_hist_tibble <- tibble::tibble(n_edges = 0,
                                  G = list(),
                                  logLG = 0)

  while(sum(G_hat) > 0){ # while the graph has edges
    if(sum(G_hat) == sum(G_complete))
      G_hist_tibble <- dplyr::bind_rows(G_hist_tibble,
                                        tibble::tibble(n_edges = sum(G_hat)/2,
                                                       G = list(G_hat),
                                                       logLG = logLG_hat))
    tmp <- step_down(sample_df, G_hat, card_A, lambda)

    if(tmp$logLG > logLG_hat) {
      G_hat <- tmp$G[[1]]
      logLG_hat <- tmp$logLG

      if(show_all_steps){
        G_hist_tibble <- dplyr::bind_rows(G_hist_tibble,
                                          tibble::tibble(n_edges = sum(G_hat)/2,
                                                         G = list(G_hat),
                                                         logLG = logLG_hat))
      }
    }else{
      break
    }
  }

  if(show_all_steps) # return all graphs selected in each step
    return(G_hist_tibble)

  tibble::tibble(G_hat = list(G_hat),
                        logLG_hat = logLG_hat)
}

