#' Add edge to an adjacency matrix
#'
#' This function adds an edge between two vertices in the given adjacency matrix.
#'
#' @param G Adjacency matrix.
#' @param v1 Vertex 1.
#' @param v2 Vertex 2.
#'
#' @return An adjacency matrix with an edge added between v1 and v2.
#' @export
#'
#' @examples
#' G_empty <- matrix(0, nrow = 5, ncol = 5)
#' add_edge(G_empty, 1, 2)
#'
add_edge <- function(G, v1, v2){
  G[v1, v2] <- 1
  G[v2, v1] <- 1
  return(G)
}


#' List possible edges to add to an adjacency matrix
#'
#' This function lists all possible edges that can be added to the given adjacency matrix to form a graph.
#'
#' @param G Adjacency matrix of a graph.
#'
#' @return A tibble with a list of all possible edges to add to graph G.
#' @export
#'
#' @examples
#' G_empty <- matrix(0, nrow = 5, ncol = 5)
#' list_possible_edges(G_empty)
#'
list_possible_edges <- function(G){
  suppressWarnings(
    G %>%
      tibble::as_tibble() %>%
      tibble::rownames_to_column("a") %>%
      tidyr::pivot_longer(names_to = "b",
                          values_to = "edge", cols = -a) %>%
      dplyr::mutate(b = as.integer(stringr::str_remove(b, "V")),
                    a = as.integer(a)) %>%
      dplyr::filter(edge == 0,
                    a < b)
  )
}


#' Step up of the forward stepwise greedy algorithm to estimate the graph
#'
#' This function implements the step up of the forward stepwise greedy algorithm to estimate the graph.
# It adds one edge to the given adjacency matrix at a time and evaluates the change in the penalized
# likelihood. It selects the edge that maximizes the increase in the likelihood and returns the updated
# adjacency matrix.
#'
#' @param sample_df Sample data.
#' @param G Adjacency matrix.
#' @param card_A Cardinality of the alphabet.
#' @param lambda Term in the penalty.
#'
#' @return A tibble with the updated adjacency matrix and the corresponding log likelihood value.
#' @export
step_up <- function(sample_df, G, card_A, lambda) {

  res <- ((list_possible_edges(G))) %>%
    dplyr::bind_cols(tibble::tibble(card_A = card_A,
                                    lambda = lambda)) %>%
    dplyr::mutate(G_plus = purrr::pmap(list(G = list(G), v1 = a, v2 = b), add_edge),
                  logLG = purrr::pmap(list(sample_data =list(sample_df),
                                           G = G_plus,
                                           card_A = card_A,
                                           lambda = lambda),
                                      get_and_sum_terms_with_penalty)) %>%
    tidyr::unnest(logLG) %>%
    dplyr::select(logLG, G_plus)

  G_plus <- (res %>%
               dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
               dplyr::pull(G_plus))[[1]]

  logLG_plus <- (res %>%
                   dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
                   dplyr::pull(logLG))[1]

  tibble::tibble(G = list(G_plus), logLG = logLG_plus)
}

#' Step up of the forward stepwise greedy algorithm to estimate the graph in parallel
#'
#' This function implements the step up of the forward stepwise greedy algorithm to estimate the graph in parallel.
# It adds one edge to the given adjacency matrix at a time and evaluates the change in the penalized
# likelihood. It selects the edge that maximizes the increase in the likelihood and returns the updated
# adjacency matrix.
#'
#' @param sample_df Sample data.
#' @param G Adjacency matrix.
#' @param card_A Cardinality of the alphabet.
#' @param lambda Term in the penalty.
#' @param clusters Number of parallel workers to use.
#'
#' @return A tibble with the updated adjacency matrix and the corresponding log likelihood value.
#' @export
step_up_parallel <- function(sample_df, G, card_A, lambda, clusters = 3) {
  future::plan(future::multisession, workers = clusters)

  possible_edges <- (list_possible_edges(G))

  result_parallel <- foreach::foreach(i = 1:nrow(possible_edges), .inorder = FALSE,
                         .options.future = list(packages=c("foreach",
                                                           "doParallel",
                                                           "magrittr",
                                                           "MixingGraph"))
                        ) %dofuture% {
    G_plus <- add_edge(G, possible_edges$a[i], possible_edges$b[i])
    logLG <- get_and_sum_terms_with_penalty(sample_data = sample_df,
                                   G = G_plus,
                                   card_A = card_A,
                                   lambda = lambda)
    tibble::tibble(logLG = logLG, G_plus = list(G_plus))
  }

  res <- tibble::enframe(result_parallel) %>%
    dplyr::select(-name) %>%
    tidyr::unnest(value)

  G_plus <- (res %>%
    dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
    dplyr::pull(G_plus))[[1]]


  logLG_plus <- (res %>%
                   dplyr::slice_max(logLG, n = 1, with_ties = FALSE) %>%
                   dplyr::pull(logLG))[1]

  future::resetWorkers(future::plan(future::sequential))

  tibble::tibble(G = list(G_plus), logLG = logLG_plus)
}

#' Implements the greedy forward stepwise selection algorithm
#'
#' This function implements the greedy forward stepwise selection algorithm to estimate the graph.
#' It iteratively adds one edge at a time to the given adjacency matrix and evaluates the change in the
#' penalized likelihood. It selects the edge that maximizes the increase in the likelihood and continues
#' until no more edges can be added without violating the acyclicity constraint of the graph.
#'
#' @param sample_df Sample data.
#' @param lambda Penalty term.
#' @param card_A Cardinality of the alphabet.
#' @param d Number of vertices in the graph.
#' @param show_all_steps Logical indicating whether all steps of the algorithm should be shown.
#'
#' @return A tibble with the estimated adjacency matrix and the corresponding penalized pseudolikelihood value.
#' @export
greedy_forward <- function(sample_df, lambda, card_A, d, show_all_steps = FALSE) {

  G_hat <- matrix(0, nrow = d, ncol = d)

  G_hist_tibble <- tibble::tibble(n_edges = 0,
                           G = list(),
                           logLG = 0)

  i = 0;

  while(sum(G_hat)/2 <= (d^2-d)/2){
    logLG_hat <- get_and_sum_terms_with_penalty(sample_df, G_hat, card_A, lambda)

    if(sum(G_hat) == 0 && show_all_steps)
      G_hist_tibble <- dplyr::bind_rows(G_hist_tibble,
                                 tibble::tibble(n_edges = sum(G_hat)/2,
                                                G = list(G_hat),
                                                logLG = logLG_hat))

    tmp <- step_up_parallel(sample_df, G_hat, card_A, lambda)

    if(tmp$logLG > logLG_hat) {
      G_hat <- tmp$G[[1]]
      logLG_hat <- tmp$logLG

      if(show_all_steps)
        G_hist_tibble <- dplyr::bind_rows(G_hist_tibble,
                                   tibble::tibble(n_edges = sum(G_hat)/2,
                                                  G = list(G_hat),
                                                  logLG = logLG_hat))
    } else{
      break
    }

    i <- i+1
    if(i > (d^2-d)/2){
      warning("Iteration reached maximum and left to avoid infinite loop")
      break
    }
  }

  if(show_all_steps) # return all graphs selected in each step
    return(G_hist_tibble)

  tibble::tibble(G_hat = list(G_hat), logLG_hat=logLG_hat)
}
