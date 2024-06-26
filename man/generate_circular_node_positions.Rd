% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{generate_circular_node_positions}
\alias{generate_circular_node_positions}
\title{Generate circular node positions}
\usage{
generate_circular_node_positions(num_nodes, radius, clockwise = TRUE)
}
\arguments{
\item{num_nodes}{Number of nodes in the circular layout.}

\item{radius}{Radius of the circle.}

\item{clockwise}{Logical indicating whether the nodes should be placed clockwise (TRUE) or counterclockwise (FALSE) on the circle.}
}
\value{
A matrix of node positions where each row represents the coordinates (x, y) of a node.
}
\description{
This function generates node positions for a circular layout with a specified number of nodes and radius.
The nodes are evenly distributed along the circumference of the circle.
}
\examples{
# Example usage:
generate_circular_node_positions(5, 10)
}
