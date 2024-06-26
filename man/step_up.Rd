% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy_fwd.R
\name{step_up}
\alias{step_up}
\title{Step up of the forward stepwise greedy algorithm to estimate the graph}
\usage{
step_up(sample_df, G, card_A, lambda)
}
\arguments{
\item{sample_df}{Sample data.}

\item{G}{Adjacency matrix.}

\item{card_A}{Cardinality of the alphabet.}

\item{lambda}{Term in the penalty.}
}
\value{
A tibble with the updated adjacency matrix and the corresponding log likelihood value.
}
\description{
This function implements the step up of the forward stepwise greedy algorithm to estimate the graph.
}
