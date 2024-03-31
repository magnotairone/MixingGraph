#' Generate Cross-Validation Data Splits
#'
#' This function generates data splits for cross-validation.
#'
#' @param original_data The original dataset to be split.
#' @param n_folds The number of folds for cross-validation.
#' @return A tibble containing the generated cross-validation data splits.
#' The returned tibble has the following columns:
#' \describe{
#'   \item{fold}{The fold number.}
#'   \item{N}{The number of samples in the fold (test data).}
#'   \item{train}{A list column containing the training data for the fold.}
#'   \item{test}{A list column containing the test data for the fold.}
#' }
#' @export
generate_cv_data <- function(original_data, n_folds) {
  res <- tidyr::tibble(fold = numeric(),
                       N = numeric(),
                       train = list(),
                       test = list())

  n <- nrow(original_data)
  folds <- (((1:n) %% n_folds) + 1) %>% sample()

  for(k in 1:n_folds){
    res <- dplyr::bind_rows(res,
                            tidyr::tibble(fold = k,
                                          N = sum(folds == k),
                                          train = list(original_data[folds != k,]),
                                          test = list(original_data[folds == k,])),)
  }
  return(res)
}

#' Generate and Save Cross-Validation Data Splits
#'
#' This function generates cross-validation data splits using the \code{generate_cv_data} function and saves the resulting data to a specified file.
#'
#' @param original_data The original dataset to be split.
#' @param n_folds The number of folds for cross-validation.
#' @param path_to_save The path where the generated data splits will be saved.
#' @param file_name The name of the file to save the data splits.
#' @param seed The random seed for reproducibility (default is 1).
#' @return None (data splits are saved to the specified file).
#' @export
generate_and_save_cv_data <- function(original_data, n_folds, path_to_save, file_name, seed = 1) {
  set.seed(seed)
  data <- generate_cv_data(original_data, n_folds)
  save(data, file = paste0(path_to_save, file_name))
}


run_algorithms_cv <- function(alg_names,
                              data, penalizing_constants, card_A,
                              path_to_partial_results, path_to_log, gamma = 1/10000) {

  if(!dir.exists(path_to_partial_results))
    dir.create(path_to_partial_results)

  folds <- unique(data$fold)

  result_cv <- tidyr::tibble(alg_name = character(),
                             c = numeric(),
                             fold = numeric(),
                             G_hat = list(),
                             logLG_pen_train = numeric(),
                             logLG_test = numeric())

  for(alg_name in alg_names){
    if (!dir.exists(paste0(path_to_partial_results, alg_name))) {
      dir.create(paste0(path_to_partial_results, alg_name))
    }

    for(penalization in penalizing_constants){
      for(k in folds) {
        df_train <- (data %>%
                       dplyr::filter(fold == k) %>%
                       dplyr::pull(train))[[1]]

        df_test <- (data %>%
                      dplyr::filter(fold == k) %>%
                      dplyr::pull(test))[[1]]

        tictoc::tic()
        if(alg_name == "fwd"){
          res <- greedy_forward(df_train,
                                lambda = penalization * log(nrow(df_train)),
                                card_A = card_A,
                                d = ncol(df_train))

        }else if(alg_name == "bcw"){
          res <- greedy_backwards(df_train,
                                  lambda = penalization * log(nrow(df_train)),
                                  card_A = card_A,
                                  d = ncol(df_train))
        }
        t2 <- tictoc::toc()

        logLG_test <- get_and_sum_terms_for_test_data(df_train,
                                                      df_test,
                                                      G_hat=res$G_hat[[1]],
                                                      card_A,
                                                      gamma)

                result_cv <- tibble(alg_name = alg_name,
                            c = penalization,
                            fold = k,
                            G_hat = res$G_hat,
                            logLG_pen_train = res$logLG_hat,
                            logLG_test = logLG_test)

        cat(paste0(alg_name, ";", penalization, ";", k, ";", t2$toc - t2$tic, ";", t2$callback_msg, ";",
                   str_split(t2$callback_msg, " ")[[1]][2], ";", res$logLG_hat, "; \n"),
            file = path_to_log, append = TRUE)

        file_path_name <- paste0(path_to_partial_results, alg_name, "/", alg_name,
                                 "_c_", penalization, "_fold_", k, ".rdata")
        save(result_cv, file = file_path_name)
      }
    }
  }
}

cv_get_data_from_folder <- function(path, folder_name){
  files <- list.files(paste0(path, folder_name))

  data <- tibble(alg_name = character(),
                 lambda = numeric(),
                 fold = numeric(),
                 G_hat = list(),
                 logLG_pen_train = numeric(),
                 logLG_test = numeric())

  for(file in files){
    load(paste0(path, folder_name, "/", file))

    splits <- (str_remove(file, ".rdata") %>%
                 str_split("_"))[[1]]

    alg_name <- splits[1]
    lambda <- splits[3] %>% as.numeric()
    fold <- splits[5] %>% as.numeric()

    data <- data %>% bind_rows(
      tibble(alg_name = alg_name,
             lambda = lambda,
             fold = fold,
             G_hat = result_cv$G_hat,
             logLG_pen_train = result_cv$logLG_pen_train,
             logLG_test = result_cv$logLG_test))
  }

  data %>% arrange(alg_name, lambda)
}

p <- function(path) {
  if(length(list.files(path)) > 1)
    stop("Theres more than one best cv result!")

  file <- list.files(path)
  load(list.files(path, full.names = T))
  res
}

# data vem do get_data_from_folder_cv
cv_get_best_lambda_from_data <- function(data) {
  data %>%
    group_by(alg_name, lambda) %>%
    summarise(logLG_cv = mean(logLG_test)) %>%
    slice_max(logLG_cv, n = 1, with_ties = TRUE) %>%
    dplyr::rename(c = lambda)
}

get_best_lambda_from_folder <- function(path, folder_name){
  cv_get_data_from_folder(path, folder_name) %>%
    filter(lambda != 0) %>%
    cv_get_best_lambda_from_data()
}

rerun_best_cv <- function(path_to_partial_results, sub_folders,
                          original_data, card_A,
                          path_to_save, path_to_log) {

  result <- tibble(path = path_to_partial_results,
                   subfolders = sub_folders) %>%
    mutate(res = pmap(list(path,  sub_folders), get_best_lambda_from_folder)) %>%
    unnest(res)

  # result
  for(i in 1:nrow(result)) {
    cat(result$alg_name[i], " c: ", result$c[i], "\n")
    run_algorithms(result$alg_name[i],
                   data = tidyr::tibble(N = nrow(original_data),
                                        data = list(original_data),
                                        sample = 1),
                   result$c[i], card_A, ncol(original_data),
                   path_to_save, path_to_log)
  }
}


generate_estimated_graph_cv <- function(path_to_folder, coords,
                                        path_to_save = NULL, height = 4, width = 6) {
  files_paths <- list.files(path_to_folder, full.names = TRUE)
  files_names <- paste0("best_cv_", list.files(path_to_folder, full.names = FALSE), ".pdf") %>%
    str_remove("_rep_1") %>%
    str_remove(".rdata")

  for (i in 1:length(files_paths)) {
    load(files_paths[i])

    g <- plot_graph_from_matrix_fixed(res$G_hat[[1]],
                                      coords, frame_border = FALSE)

    if(!is.null(path_to_save)){
      file_name <- paste0(path_to_save, files_names[i])
      cat(file_name, "\n")

      export_graphic_to_pdf(
        file_name,
        g,
        height, width
      )

    } else {
      print(g)
    }
  }
}
