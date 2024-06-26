% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudolikelihood.R
\name{update_sum_terms}
\alias{update_sum_terms}
\title{Update sum terms}
\usage{
update_sum_terms(terms, sample_data, G, vertices)
}
\arguments{
\item{terms}{Result of *get_sum_terms* function.}

\item{sample_data}{Data sample.}

\item{G}{Adjacency matrix of the given graph.}

\item{vertices}{Vertices to update their terms in the sum.}
}
\value{
Updated sum terms.
}
\description{
This function updates the sum terms for the given vertices in the pseudolikelihood calculation.
It takes the adjacency matrix of the graph, vertices to update, the result of the *get_sum_terms* function,
and the sample data as input, and returns the updated sum terms.
}
