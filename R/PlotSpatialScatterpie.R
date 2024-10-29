#' @title Draw the scatterpie plot for spatial transcriptomics data
#'
#' @description This function can draw the scatterpie plot for spatial transcriptomics data.
#'
#' @param st_vis spatial transcriptomic data: Seurat object
#' @param cell_types_all a vecter or list data includes the cell type information
#' @param img_path the corresponding image path for st_vis
#' @param cols optional: set the color for each cell type in the pie plot. By default is 'NULL'.
#' @param cell_types_interest optional: Used cell types for drawing the pie plot. By default is 'NULL', which means use all cell types.
#' @param slice optional: set the slice used for drawing the pie plot. By default is 'NULL'.
#' @param scatterpie_alpha optional: set the diaphaneity/alpha of the pie. By default is '1'.
#' @param pie_scale optional: set the size of the pie. By default is '1'.
#'
#' @importFrom grid rasterGrob
#' @importFrom grid unit
#' @importFrom png readPNG
#' @importFrom jpeg readJPEG
#' @importFrom dplyr %>%
#'
#' @return a ggplot2 object
#' @export


spatial_scatterpie <- function (st_vis, cell_types_all, img_path,cols = NULL, cell_types_interest = NULL,
                                slice = NULL, scatterpie_alpha = 1, pie_scale = 1)
{
  if (!is(st_vis, "Seurat"))
    stop("ERROR: st_vis must be a Seurat object!")
  if (!is(cell_types_all, "vector"))
    stop("ERROR: cell_types_all must be a vector/list object!")
  if (!is.character(img_path))
    stop("ERROR: must be a character string!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest)))
    stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is.numeric(scatterpie_alpha))
    stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale))
    stop("ERROR: pie_scale must be numeric between 0 and 1!")
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  # suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))
  suppressMessages(require(png))
  suppressMessages(require(jpeg))
  suppressMessages(require(grid))
  metadata_ds <- data.frame(st_vis@meta.data)
  colnames(metadata_ds) <- colnames(st_vis@meta.data)
  if (is.null(cell_types_interest)) {
    cell_types_interest <- cell_types_all
  }
  if (!all(cell_types_all %in% cell_types_interest)) {
    metadata_ds <- metadata_ds %>% tibble::rownames_to_column("barcodeID") %>%
      dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest,
                                           drop = FALSE])) %>% dplyr::filter(rsum != 0) %>%
      dplyr::select("barcodeID") %>% dplyr::left_join(metadata_ds %>%
                                                        tibble::rownames_to_column("barcodeID"), by = "barcodeID") %>%
      tibble::column_to_rownames("barcodeID")
  }
  if (is.null(slice) | (!is.null(slice) && !slice %in% names(st_vis@images))) {
    slice <- names(st_vis@images)[1]
    warning(sprintf("Using slice %s", slice))
  }
  spatial_coord <- data.frame(st_vis@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow*st_vis@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol*st_vis@images[[slice]]@scale.factors$lowres) %>%
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"), by = "barcodeID")

  img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
  if (img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  }
  else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }
  img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  if (is.null(cols))
  {
    scatterpie_plt <-
      suppressMessages(
        ggplot2::ggplot() + ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)
        ) +
          scatterpie::geom_scatterpie(
            data = spatial_coord,
            ggplot2::aes(x = imagecol_scaled,
                         y = imagerow_scaled),
            cols = cell_types_all,
            color = NA,
            alpha = scatterpie_alpha,
            pie_scale = pie_scale
          ) +
          ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img),
                                                     0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11,
                                                                                                                 rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(
                                                                                                                   ratio = 1,
                                                                                                                   xlim = NULL,
                                                                                                                   ylim = NULL,
                                                                                                                   expand = TRUE,
                                                                                                                   clip = "on"
                                                                                                                 )
      )
  }
  else
  {
    scatterpie_plt <- suppressMessages(ggplot2::ggplot() +
                                         ggplot2::annotation_custom(grob = img_grob,xmin = 0,
                                                                    xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                                         scatterpie::geom_scatterpie(data = spatial_coord,
                                                                     mapping = ggplot2::aes(x = imagecol_scaled,y = imagerow_scaled),
                                                                     cols = cell_types_all, color = NA,
                                                                     alpha = scatterpie_alpha, pie_scale = pie_scale) +
                                         scale_fill_manual(values=cols)+
                                         ggplot2::scale_y_reverse() +
                                         ggplot2::ylim(nrow(img), 0) +
                                         ggplot2::xlim(0, ncol(img)) +
                                         cowplot::theme_half_open(11, rel_small = 1) +
                                         ggplot2::theme_void() +
                                         ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))

  }
  return(scatterpie_plt)
}
