#' Plot Dimensionality Reduction
#'
#' This function performs dimensionality reduction on gene expression data
#' and visualizes the results using `plotly`, integrating metadata information
#' whenever possible.
#' It supports multiple methods, including PCA, MDS, t-SNE, and UMAP.
#'
#' @param compendium A data frame or matrix where rows represent genes
#' and columns represent samples.
#' @param metadata A data frame containing metadata for the samples,
#' where each row corresponds to a sample in `compendium`,
#' if metadata for that sample are available.
#' If NULL tooltips will just display the ID of the sample contained in 
#' `compendium`
#' @param color_by A string specifying the column name in `metadata`
#' to be used for coloring points in the plot. MUST be NULL if `metadata` is 
#' NULL
#' @param method A string specifying the dimensionality reduction method.
#' Options are "PCA", "MDS", "t-SNE", or "UMAP". Defaults to "PCA".
#'
#' @return A `plotly` object visualizing the dimensionality reduction result,
#' with interactive tooltips, showing metadata information.
#'
#' @details
#' - PCA (Principal Component Analysis): Computes principal components and 
#' explained variance.
#' - MDS (Multidimensional Scaling): Uses pairwise distances between samples for
#'  dimensionality reduction.
#' - t-SNE (t-Distributed Stochastic Neighbor Embedding): Optimizes a similarity
#'  measure for high-dimensional data visualization.
#' - UMAP (Uniform Manifold Approximation and Projection): Aims to preserve 
#' local and global structure of data.
#' @importFrom magrittr %>%
#' @importFrom stats prcomp dist cmdscale setNames
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr mutate left_join
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom plotly plot_ly layout
#' @export
plot_dimensionality_reduction <- function(compendium, metadata = NULL, color_by = NULL, method = "PCA"){

  if(is.null(metadata)&&!is.null(color_by)){
    stop("When metadata is NULL, color_by MUST also be NULL")
  }

  `%>%` <- magrittr::`%>%`
  data_transposed <- t(compendium)

  reduction_result <- NULL
  explained_variance <- NULL
  axis_labels <- c("Dimension 1", "Dimension 2")

  if (method == "PCA") {
    pca <- stats::prcomp(data_transposed, center = TRUE, scale. = TRUE)
    reduction_result <- suppressWarnings(tibble::as_tibble(pca$x[, 1:2]))
    explained_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
    colnames(reduction_result) <- c("PC1", "PC2")
    axis_labels <- c(paste0("PC1 (", explained_variance[1], "%)"), paste0("PC2 (", explained_variance[2], "%)"))
  } else if (method == "MDS") {
    distance_matrix <- stats::dist(data_transposed)
    mds <- stats::cmdscale(distance_matrix, k = 2)
    reduction_result <- suppressWarnings(tibble::as_tibble(mds))
    colnames(reduction_result) <- c("MDS1", "MDS2")
    axis_labels <- c("MDS1", "MDS2")
  } else if (method == "t-SNE") {
    tsne <- Rtsne::Rtsne(data_transposed, dims = 2, perplexity = 30, verbose = FALSE, max_iter = 500)
    reduction_result <- suppressWarnings(tibble::as_tibble(tsne$Y))
    colnames(reduction_result) <- c("t-SNE1", "t-SNE2")
    axis_labels <- c("t-SNE1", "t-SNE2")
  } else if (method == "UMAP") {
    umap_result <- umap::umap(data_transposed, n_neighbors = 15, n_components = 2, metric = "euclidean")
    reduction_result <- suppressWarnings(tibble::as_tibble(umap_result$layout))
    colnames(reduction_result) <- c("UMAP1", "UMAP2")
    axis_labels <- c("UMAP1", "UMAP2")
  } else {
    stop("Method not supported. Please choose among 'PCA', 'MDS', 't-SNE', or 'UMAP'.")
  }

  reduction_result <- reduction_result %>% dplyr::mutate(Sample = rownames(data_transposed))

  if (is.null(metadata) && is.null(color_by)) {
    reduction_result <- reduction_result %>% dplyr::mutate(Tooltip = paste0("run_accession: ",colnames(compendium)))

    fig <- plotly::plot_ly(
      data = reduction_result,
      x = ~ .data[[colnames(reduction_result)[1]]],
      y = ~ .data[[colnames(reduction_result)[2]]],
      type = 'scatter', mode = 'markers',
      marker = list(size = 10, color = '#636EFA'),
      text = ~Tooltip,
      hoverinfo = "text"
    ) %>% plotly::layout(
      title = paste(method, 'on Conditions'),
      xaxis = list(title = axis_labels[1]),
      yaxis = list(title = axis_labels[2]),
      legend = list(orientation = 'h', x = 0, y = -0.2)
    )

    fig <- htmlwidgets::onRender(fig, "
  function(el, x) {
    el.on('plotly_selected', function(data) {
      if (data && data.points) {
        var selected_accessions = data.points.map(pt => {
          var tooltip = pt.text || '';
          var accession_match = tooltip.match(/run_accession: ([^<]+?)(?=<br>|$)/);
          return accession_match ? accession_match[1] : null;
        }).filter(accession => accession !== null);
        Shiny.onInputChange('selected_run_accession', selected_accessions);
      } else {
        Shiny.onInputChange('selected_run_accession', null);
      }
    });
  }
");

    return(fig)
  }

  if (!"Sample" %in% colnames(metadata)) {
    metadata <- metadata %>% tibble::rownames_to_column(var = "Sample")
  }

  if (!color_by %in% colnames(metadata) && !is.null(color_by)) {
    stop(paste("Column", color_by, "does not exist in 'metadata' data.frame"))
  }

  if (is.null(color_by)) {
    tooltip_text <- apply(metadata, 1, function(row) {
      to_keep <- c("study_accession", "tax_id", "scientific_name", "instrument_model", "sample_title")
      paste(
        paste0("run_accession: ", row["run_accession"]),
        paste(names(row[to_keep]), row[to_keep], sep = ": ", collapse = "<br>"),
        sep = "<br>"
      )
    })

    reduction_data <- reduction_result %>% dplyr::left_join(metadata, by = "Sample")

    fig <- plotly::plot_ly(
      data = reduction_data,
      x = ~ .data[[colnames(reduction_result)[1]]],
      y = ~ .data[[colnames(reduction_result)[2]]],
      type = 'scatter', mode = 'markers',
      marker = list(size = 10, color = '#636EFA'),
      text = ~tooltip_text[match(Sample, metadata$Sample)],
      hoverinfo = "text"
    ) %>% plotly::layout(
      title = paste(method, 'on Conditions'),
      xaxis = list(title = axis_labels[1]),
      yaxis = list(title = axis_labels[2]),
      legend = list(orientation = 'h', x = 0, y = -0.2)
    )

    fig <- htmlwidgets::onRender(fig, "
    function(el, x) {
      el.on('plotly_selected', function(data) {
        if (data && data.points) {
          var selected_accessions = data.points.map(pt => {
            var tooltip = pt.text || '';
            var accession_match = tooltip.match(/run_accession: ([^<]+?)(?=<br>|$)/);
            return accession_match ? accession_match[1] : null;
          }).filter(accession => accession !== null);
          Shiny.onInputChange('selected_run_accession', selected_accessions);
        } else {
          Shiny.onInputChange('selected_run_accession', null);
        }
      });
    }
    ");

    return(fig)
  }

  reduction_data <- reduction_result %>%
    dplyr::left_join(metadata, by = "Sample")

  unique_groups <- unique(reduction_data[[color_by]])
  num_groups <- length(unique_groups)

  palette_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, num_groups), "Set1"))(num_groups)
  color_mapping <- stats::setNames(palette_colors, unique_groups)

  tooltip_text <- apply(metadata, 1, function(row) {
    to_keep <- c("study_accession", "tax_id", "scientific_name", "instrument_model", "sample_title")
    paste(
      paste0("run_accession: ", row["run_accession"]),
      paste(names(row[to_keep]), row[to_keep], sep = ": ", collapse = "<br>"),
      sep = "<br>"
    )
  })

  fig <- plotly::plot_ly(
    data = reduction_data,
    x = ~ .data[[colnames(reduction_result)[1]]],
    y = ~ .data[[colnames(reduction_result)[2]]],
    type = 'scatter', mode = 'markers',
    color = ~.data[[color_by]],
    colors = palette_colors,
    text = ~tooltip_text[match(Sample, metadata$Sample)],
    marker = list(size = 10),
    hoverinfo = "text"
  ) %>% plotly::layout(
    title = paste(method, 'on Conditions'),
    xaxis = list(title = axis_labels[1]),
    yaxis = list(title = axis_labels[2]),
    legend = list(
      title = list(text = color_by),
      orientation = 'h',
      x = 0, y = -0.2,
      font = list(size = 10),
      itemwidth = 30,
      tracegroupgap = 5
    )
  )

  fig <- htmlwidgets::onRender(fig, "
  function(el, x) {
    el.on('plotly_selected', function(data) {
      if (data && data.points) {
        var selected_accessions = data.points.map(pt => {
          var tooltip = pt.text || '';
          var accession_match = tooltip.match(/run_accession: ([^<]+?)(?=<br>|$)/);
          return accession_match ? accession_match[1] : null;
        }).filter(accession => accession !== null);
        Shiny.onInputChange('selected_run_accession', selected_accessions);
      } else {
        Shiny.onInputChange('selected_run_accession', null);
      }
    });
  }
");

  fig <- fig %>% layout(
    legend = list(
      x = -0.3,
      y = -0.7,
      xanchor = "down"
    ),
    margin = list(l = 50, r = 150, b = 50, t = 50)
  )

  return(fig)
}
