#' Randomly choose a neighbor for a given graph.
#'
#' This function randomly selects a neighbor for a given graph by flipping the state of one randomly chosen edge.
#'
#' @param G The graph from which to choose a neighbor.
#' @param d The number of nodes in G. Defaults to the number of rows in G.
#'
#' @return The adjacency matrix of the graph with a randomly chosen neighbor.
#' @export
#'
#' @examples
#' G <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
#' set.seed(1)
#' choose_neighbor(G)
choose_neighbor <- function(G, d = nrow(G)){

  V <- sample(d, 2)

  new_state <- as.numeric(!G[V[1],V[2]])

  G[V[1],V[2]] <- new_state
  G[V[2],V[1]] <- new_state

  G
}

#' Identify differences in the edges of graphs G1 and G2.
#'
#' This function identifies the differing edge between two graphs G1 and G2. It assumes that G1 and G2 are neighbors, meaning they differ by only one edge.
#'
#' @param G1 The adjacency matrix of the first graph.
#' @param G2 The adjacency matrix of the second graph.
#'
#' @return A numeric vector representing the position of the differing edge.
#' @export
#'
#' @examples
#' G1 <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
#' set.seed(1)
#' G2 <- choose_neighbor(G1)
#' identify_diff(G1, G2)
identify_diff <- function(G1, G2){
  which(G1 != G2, arr.ind = TRUE)[1,]
}

#' Simulated annealing algorithm
#'
#' This function implements the simulated annealing algorithm for estimating the adjacency matrix of a graph.
#'
#' @param sample_data Sample data used for estimation.
#' @param lambda Penalty term used in the estimation.
#' @param card_A Cardinality of the alphabet.
#' @param G_0 Initial guess for the adjacency matrix.
#' @param N Number of iterations of the simulated annealing algorithm.
#' @param C Constant parameter in the simulated annealing algorithm.
#' @param include_logLG If TRUE, returns both the estimated adjacency matrix and the corresponding penalized pseudolikelihood value.
#'
#' @return If `include_logLG` is FALSE, returns the estimated adjacency matrix. If TRUE, returns a tibble containing both the estimated adjacency matrix (`G_hat`) and the corresponding penalized pseudolikelihood (`logLG_hat`).
#' @export
#'
SimAn <- function(sample_data, lambda, card_A, G_0, N, C,
                  include_logLG = FALSE){
  G_current <- G_0
  terms_list <- get_sum_terms(sample_data, G_current)
  logLG_hat <- sum_terms(terms_list) - penalty(lambda, card_A, G_current)

  G_hat <- G_current
  logLG_current <- logLG_hat

  for(i in 1:N){

    G <- choose_neighbor(G_current)
    W <- identify_diff(G_current, G)

    terms_list <- update_sum_terms(terms_list, sample_data, G, W)
    logLG_tmp <- sum_terms(terms_list) - penalty(lambda, card_A, G)

    eta <- C * log(1+i)

    if(logLG_tmp > logLG_current) {
      G_current <- G
      logLG_current <- logLG_tmp

      if(logLG_tmp > logLG_hat) {
        G_hat <- G_current
        logLG_hat <- logLG_current
      }
    }else {
      p <- exp(eta * (logLG_tmp - logLG_current))
      u <- runif(1)

      if(u < p) {
        G_current <- G
        logLG_current <- logLG_tmp
      }
    }
  }

  if(include_logLG) {
    return(dplyr::tibble(G_hat = list(G_hat), logLG_hat = logLG_hat))
  }else{
    return(G_hat)
  }
}
