#' Get sum terms of the pseudolikelihood
#'
#' This function calculates the sum terms of the pseudolikelihood for each vertex in the graph.
#' It computes the term for each vertex based on the given sample data and adjacency matrix.
#'
#' @param sample_df Data sample.
#' @param G Adjacency matrix of the given graph.
#'
#' @return A list of the terms for each vertex in the graph. Each item is a tibble.
#' @export
get_sum_terms <- function(sample_df, G) {
  result <- list()

  for (v in 1:ncol(G)) {
    if (length(which(G[v, ] == 1)) == 0) {# W is empty
      V <- paste("V", v, sep = "")

      a_v_W <- sample_df %>%
        dplyr::group_by_at(c(V)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(N_a_w = nrow(sample_df),
               pi = N_a_v_W * log(N_a_v_W / N_a_w))

      result[[v]] <- a_v_W
      next
    }

    W <- paste("V", which(G[v, ] == 1), sep = "")
    V <- paste("V", v, sep = "")

    a_v_W <- sample_df %>%
      dplyr::group_by_at(c(V, W)) %>%
      dplyr::count(name = "N_a_v_W") %>%
      dplyr::ungroup()

    suppressMessages(
      suppressWarnings(
        a_v <- a_v_W %>%
          dplyr::group_by_at(W) %>%
          dplyr::summarize(N_a_w = sum(N_a_v_W)) %>%
          dplyr::ungroup()
      ))

    suppressMessages(
      suppressWarnings(
        result[[v]] <- a_v_W %>%
          dplyr::left_join(a_v) %>%
          dplyr::mutate(pi = N_a_v_W * log(N_a_v_W / N_a_w))
      ))
  }

  result
}

get_sum_terms_parallel <- function(sample_df, G) {
  result <- list()

  future::plan(future::multisession, workers = 3)
  x <- foreach::foreach(v = 1:ncol(G),
                        .inorder = FALSE,
                        .options.future = list(packages=c("foreach", "doParallel",
                                                          "magrittr", "magnoPhD"))
                        ) %dofuture% {
    if (length(which(G[v, ] == 1)) == 0) {# W is empty
      V <- paste("V", v, sep = "")

      a_v_W <- sample_df %>%
        dplyr::group_by_at(c(V)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(N_a_w = nrow(sample_df),
                      pi = N_a_v_W * log(N_a_v_W / N_a_w))

      result[[v]] <- a_v_W
    }else{
      W <- paste("V", which(G[v, ] == 1), sep = "")
      V <- paste("V", v, sep = "")

      a_v_W <- sample_df %>%
        dplyr::group_by_at(c(V, W)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup()

      suppressMessages(
        suppressWarnings(
          a_v <- a_v_W %>%
            dplyr::group_by_at(W) %>%
            dplyr::summarize(N_a_W = sum(N_a_v_W)) %>%
            dplyr::ungroup()
        ))

      suppressMessages(
        suppressWarnings(
          result[[v]] <- a_v_W %>%
            dplyr::left_join(a_v) %>%
            dplyr::mutate(pi = N_a_v_W * log(N_a_v_W / N_a_w))
        ))
    }
  }
  future::resetWorkers(future::plan(future::sequential))

  return(result)
}

#' sum_terms
#'
#' @param terms result of *get_sum_terms* function
#'
#' @return a numeric value representing the sum of the terms
#' @export
#'
sum_terms <- function(terms) {
  purrr::map_dbl(terms, function(x)
    sum(x$pi)) %>%
    sum()
}

#' Update sum terms
#'
#' This function updates the sum terms for the given vertices in the pseudolikelihood calculation.
#' It takes the adjacency matrix of the graph, vertices to update, the result of the *get_sum_terms* function,
#' and the sample data as input, and returns the updated sum terms.
#'
#' @param G Adjacency matrix of the given graph.
#' @param vertices Vertices to update their terms in the sum.
#' @param terms Result of *get_sum_terms* function.
#' @param sample_data Data sample.
#'
#' @return Updated sum terms.
#' @export
update_sum_terms <- function(terms, sample_data, G, vertices) {
  for (v in vertices) {
    if (length(which(G[v, ] == 1)) == 0) {# if W is empty
      V <- paste("V", v, sep = "")

      a_v_W <- sample_data %>%
        dplyr::group_by_at(c(V)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(N_a_W = nrow(sample_data),
               pi = N_a_v_W * log(N_a_v_W / N_a_W))

      terms[[v]] <- a_v_W
      next

    }

    W <- paste("V", which(G[v, ] == 1), sep = "")
    V <- paste("V", v, sep = "")

    a_v_W <- sample_data %>%
      dplyr::group_by_at(c(V, W)) %>%
      dplyr::count(name = "N_a_v_W") %>%
      dplyr::ungroup()

    suppressMessages(
      suppressWarnings(
        a_v <- a_v_W %>%
          dplyr::group_by_at(W) %>%
          dplyr::summarize(N_a_W = sum(N_a_v_W)) %>%
          dplyr::ungroup()
      ))

    suppressMessages(
      suppressWarnings(
        terms[[v]] <- a_v_W %>%
          dplyr::left_join(a_v) %>%
          dplyr::mutate(pi = N_a_v_W * log(N_a_v_W / N_a_W))
      ))
  }
  terms
}

#' Penalty for the likelihood
#'
#' This function calculates the penalty term for the likelihood. The penalty considered is lambda * |A|^|G|,
#' where lambda is the penalty coefficient, |A| is the cardinality of the alphabet, and |G| is the number of
#' edges in the graph represented by the adjacency matrix G.
#'
#' @param lambda Penalty coefficient.
#' @param card_A Cardinality of the alphabet.
#' @param G Adjacency matrix of the given graph.
#'
#' @return Penalty term for the likelihood.
#' @export
penalty <- function(lambda, card_A, G) {
  lambda * sum(card_A ^ colSums(G))
}

#' Get and sum terms of the pseudolikelihood function
#'
#' This function calculates and sums the terms of the pseudolikelihood function for a given graph represented
#' by the adjacency matrix G and the sample data. It first obtains the terms using the *get_sum_terms* function
#' and then sums them using the *sum_terms* function.
#'
#' @param sample_data Sample data.
#' @param G Adjacency matrix of the graph.
#'
#' @return Sum of the terms of the pseudolikelihood function.
get_and_sum_terms <- function(sample_data, G){
  get_sum_terms(sample_data, G) %>%
    sum_terms()
}

#' Get and sum terms of the penalized pseudolikelihood function
#'
#' This function calculates and sums the terms of the penalized pseudolikelihood function for a given graph
#' represented by the adjacency matrix G, the sample data, the cardinality of the alphabet, and the lambda term
#' in the penalty. It first obtains the terms using the *get_sum_terms* function, then calculates the penalty
#' using the *penalty* function, and finally subtracts the penalty from the sum of the terms.
#'
#' @param sample_data Sample data.
#' @param G Adjacency matrix of the graph.
#' @param card_A Cardinality of the alphabet.
#' @param lambda Lambda term in the penalty.
#'
#' @return Sum of the terms of the penalized pseudolikelihood function.
#' @export
get_and_sum_terms_with_penalty <- function(sample_data, G, card_A, lambda){
  terms <- get_sum_terms(sample_data, G)
  return(sum_terms(terms) - penalty(lambda, card_A, G))
}

get_and_sum_terms_with_penalty_parallel <- function(sample_data, G, card_A, lambda){
  (get_sum_terms_parallel(sample_data, G) %>%
     sum_terms()) - penalty(lambda, card_A, G)
}

get_estimated_log_probs <- function(sample_df, G_hat, card_A, gamma = 1/10000){
  result <- list()

  for (v in 1:ncol(G_hat)) {
    if (length(which(G_hat[v, ] == 1)) == 0) {# W is empty
      V <- paste("V", v, sep = "")
      W <- character()
      a_v_W <- sample_df %>%
        dplyr::group_by_at(c(V)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(N_a_W = nrow(sample_df),
                      pi_hat = (N_a_v_W / N_a_W))

      result[[v]] <- a_v_W %>% dplyr::select(V, pi_hat)

    }else{
      W <- paste("V", which(G_hat[v, ] == 1), sep = "")
      V <- paste("V", v, sep = "")

      a_v_W <- sample_df %>%
        dplyr::group_by_at(c(V, W)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup()

      suppressMessages(
        suppressWarnings(
          a_W <- a_v_W %>%
            dplyr::group_by_at(W) %>%
            dplyr::summarize(N_a_W = sum(N_a_v_W)) %>%
            dplyr::ungroup()
        ))

      suppressMessages(
        suppressWarnings(
          result[[v]] <- a_v_W %>%
            dplyr::left_join(a_W) %>%
            dplyr::mutate(pi_hat = (N_a_v_W / N_a_W)) %>%
            dplyr::select(all_of(c(V, W)), pi_hat)
        ))
    }

    # check if there are unobserved words
    n_observed_configs <- nrow(result[[v]])
    n_possible_configs <- card_A^(ncol(result[[v]])-1)

    # add word that has not been observed
    if(n_observed_configs < n_possible_configs) {
      n_nodes <- length(c(V,W))
      m <- matrix(rep(0:(card_A-1), n_nodes), ncol = n_nodes, byrow = F)
      colnames(m) <- c(V,W)

      suppressMessages(
        suppressWarnings(
          missing_combination <- expand.grid(m %>% tidyr::as_tibble()) %>%
            tidyr::as_tibble() %>%
            dplyr::anti_join(result[[v]]) %>%
            dplyr::mutate(pi_hat = (gamma),
                          added = TRUE)
        ))

      result[[v]] <-
        result[[v]] %>%
        mutate(added = FALSE) %>%
        dplyr::bind_rows(missing_combination) %>%
        group_by_at(W) %>%
        mutate(n = n(),
               observed = n()-sum(added),
               non_observed = sum(added)) %>%
        mutate(pi_hat = if_else(!added,
                                 pi_hat - gamma * non_observed / observed,
                                 pi_hat)) %>%
        mutate(pi_hat = ifelse(sum(added)==n, 1/n, pi_hat)) %>%
        arrange_at(W)
    }

    result[[v]] <- result[[v]] %>%
        mutate(log_pi = log(pi_hat)) %>%
      mutate(pi = exp(log_pi))
  }

  result
}

get_and_sum_terms_for_test_data <- function(df_train, df_test, G_hat, card_A, gamma = 1/10000) {
  tidyr::tibble(df = sum_terms_for_test_data(df_train, df_test, G_hat, card_A, gamma)) %>%
    dplyr::mutate(sum = purrr::map_dbl(df, ~ sum(.$pi))) %>%
    dplyr::summarise(sum = sum(sum)) %>%
    dplyr::pull(sum)
}

sum_terms_for_test_data <- function(df_train, df_test, G_hat, card_A, gamma = 1/10000){
  estimated_probs_from_train <- get_estimated_log_probs(df_train, G_hat, card_A, gamma)
  test_data_counts <- list()
  result <- list()

  for (v in 1:ncol(G_hat)) {
    if (length(which(G_hat[v, ] == 1)) == 0) {# if W is empty
      V <- paste("V", v, sep = "")
      W <- " "

      a_v_W <- df_test %>%
        dplyr::group_by_at(c(V)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(N_a_W = nrow(df_test))

      test_data_counts[[v]] <- a_v_W

    }else{
      W <- paste("V", which(G_hat[v, ] == 1), sep = "")
      V <- paste("V", v, sep = "")

      a_v_W <- df_test %>%
        dplyr::group_by_at(c(V, W)) %>%
        dplyr::count(name = "N_a_v_W") %>%
        dplyr::ungroup()

      suppressMessages(
        suppressWarnings(
          a_W <- a_v_W %>%
            dplyr::group_by_at(W) %>%
            dplyr::summarize(N_a_W = sum(N_a_v_W)) %>%
            dplyr::ungroup()
        ))

      suppressMessages(
        suppressWarnings(
          test_data_counts[[v]] <- a_v_W %>%
            dplyr::left_join(a_W)
        ))
    }

    suppressMessages(
      suppressWarnings(
        result[[v]] <- test_data_counts[[v]] %>%
          dplyr::left_join(estimated_probs_from_train[[v]]) %>%
          dplyr::mutate(pi = N_a_v_W * log_pi) %>%
          dplyr::select(dplyr::any_of(c(V, W, "pi", "N_a_v_W")))
      ))
  }

  result
}
