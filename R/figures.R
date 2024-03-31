#' Generate circular node positions
#'
#' This function generates node positions for a circular layout with a specified number of nodes and radius.
#' The nodes are evenly distributed along the circumference of the circle.
#'
#' @param num_nodes Number of nodes in the circular layout.
#' @param radius Radius of the circle.
#' @param clockwise Logical indicating whether the nodes should be placed clockwise (TRUE) or counterclockwise (FALSE) on the circle.
#'
#' @return A matrix of node positions where each row represents the coordinates (x, y) of a node.
#' @export
#'
#' @examples
#' # Example usage:
#' generate_circular_node_positions(5, 10)
generate_circular_node_positions <- function(num_nodes, radius, clockwise = TRUE) {
  angles <- seq(0, 2*pi, length.out = num_nodes+1)[-1]  # Start from 0 for clockwise order, or 2*pi for counterclockwise order

  if (clockwise)
    angles <- rev(angles)

  x <- radius * cos(angles)
  y <- radius * sin(angles)

  node_positions <- cbind(x, y)

  return(node_positions)
}


#' Plot Graph from Matrix with Fixed Layout
#'
#' This function plots a graph from an adjacency matrix and node coordinates, using a fixed layout.
#'
#' @param G adjacency matrix that represents the graph G
#' @param node_size node size in the plot
#' @param title plot title
#' @param edge_width edge width in the plot
#' @param edge_color edge color in the plot
#' @param coords matrix of coordinates of the vertices in the plot
#' @param frame_border should frame border be shown?
#'
#' @importFrom extrafont loadfonts
#'
#' @return the graph created from the adjacency matrix
#' @export
#'
#' @examples
#' coords <-  matrix(c(1,2, 3,1, 5,2, 4,3, 2,3), ncol = 2, byrow = TRUE)
#' G <- matrix(1, nrow = 5, ncol = 5) - diag(1, 5)
#' # plot_graph_from_matrix(G, title = "Fixed random graph", coords = coords)
plot_graph_from_matrix <-
  function(G,
           coords,
           node_size = 4,
           title = NULL,
           edge_width = 2,
           edge_color = "darkgray",
           frame_border = TRUE,
           node_font_size = 5,
           label_padding = 0.5,
           alpha_label = 1) {
    if (sum(G) == 0) { # empty
      G <- matrix(1, nrow = nrow(G), ncol = ncol(G)) -
        diag(1, nrow = nrow(G), ncol = ncol(G))
      edge_color = "white"
    }

    if(is.null(colnames(G))){

      graph <- as.data.frame(G) %>%
        dplyr::mutate(a = 1:nrow(G)) %>%
        tidyr::pivot_longer(
          names_to = "b",
          values_to =
            "edge",
          cols = -a
        ) %>%
        dplyr::mutate(b = as.integer(stringr::str_sub(b, 2))) %>%
        dplyr::filter(edge == 1, a <
                        b) %>%
        rbind(data.frame(
          a = 1:nrow(G),
          b = 1:nrow(G),
          edge = 1
        )) %>%
        dplyr::arrange(a, b) %>%
        igraph::graph_from_data_frame()

    }else{
      graph <- tibble::as_tibble(G, rownames = "a") %>%
        tidyr::pivot_longer(
          names_to = "b",
          values_to =
            "edge",
          cols = -a
        ) %>%
        dplyr::filter(edge == 1) %>%
        dplyr::arrange(a, b) %>%
        igraph::graph_from_data_frame(vertices = tibble::tibble(name = colnames(G)))

    }

    if (frame_border) {
      borda <- ggplot2::element_rect(
        fill = NA,
        color = "black",
        size = 0.5,
        linetype = "solid"
      )
    } else{
      borda <- ggplot2::element_blank()
    }

    if (!is.matrix(coords)) {
      return(errorCondition("Must input the vertices coordinates"))
    }

    lay <-
      ggraph::create_layout(graph,
                            layout = "manual",
                            x = coords[, 1],
                            y = coords[, 2])

    g <- lay %>%
      ggraph::ggraph() +
      ggraph::geom_edge_link(colour = edge_color, edge_width = edge_width) +
      ggraph::geom_node_label(
        ggplot2::aes(label = name),
        repel = F ,
        size = node_size,
        label.padding = ggplot2::unit(label_padding, "lines"),
        family = "lmroman",
        alpha = alpha_label
      ) +
      ggraph::theme_graph() + #base_family = "lmroman") +
      ggplot2::scale_y_continuous(expand = c(0.2, 0.2)) +
      ggplot2::scale_x_continuous(expand = c(0.2, 0.2))

    if (!is.null(title)) {
      g <- g +
        ggplot2::ggtitle(title)
    }

    g +
      ltxplot::theme_latex() +
      ggraph::theme_graph() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(rep(.5,4)), # margin(10, 10, 10, 10)
        strip.text = ggplot2::element_text(family = "lmroman"),
        plot.caption = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.2)),
        panel.border = borda,
        plot.title = ggplot2::element_text(size = 10),
      )
  }
