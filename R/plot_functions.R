#' Plot prototype EEMs of a self-organising map
#'
#' @description The prototype EEMs (codebook vectors) of all units or selected units are plotted. If all units are plotted, the facets are arranged like the units on the map (top row of the map first).
#'
#' @param data object of class eemsom, see \code{\link{eemsom}}
#' @param fill_max sets the maximum fluorescence value for the colour scale, see \code{\link[staRdom]{ggeem}}
#' @param units optional numeric vector of units to be plotted
#' @param denormalise logical, see \code{\link{eemsom_codes}}
#' @param ... parameters passed on to \code{\link[staRdom]{ggeem}}, e.g. \code{colpal} or \code{contour}
#'
#' @details The facet titles show the unit and the number of samples mapped to this unit. On a hexagonal map, every second row is shifted by half a unit, which cannot be shown by the facets.
#'
#' Ejarque-Gonzalez and Butturini (2014) used the prototype EEMs to describe the characteristic fluorescence of the regions of the map.
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @references Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @seealso \code{\link[staRdom]{ggeem}}, \code{\link{eemsom_codes}}
#' @family plots
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100, seed = 42)
#' ggeem(som)
#' ggeem(som, units = c(1, 2), contour = TRUE)
ggeem.eemsom <- function(data, fill_max = FALSE, units = NULL, denormalise = FALSE, ...){
  codes <- eemsom_codes(data, denormalise = denormalise)
  pts <- data$som$grid$pts
  n <- tabulate(data$som$unit.classif, nbins = nrow(pts))
  labs <- paste0(eem_names(codes), " (n = ", n, ")")
  all_units <- is.null(units)
  if(all_units) units <- seq_len(nrow(pts))
  if(!all(units %in% seq_len(nrow(pts)))) stop("units need to be numbers between 1 and ", nrow(pts), "!")
  ord <- order(-pts[, 2], pts[, 1])
  table <- lapply(units, function(u){
    as.data.frame(codes[[u]]) %>%
      mutate(sample = labs[u])
  }) %>%
    bind_rows() %>%
    mutate(sample = factor(sample, levels = labs[ord[ord %in% units]]))
  plot <- ggeem(table, fill_max = fill_max, ...)
  if(all_units) plot <- plot + facet_wrap(~ sample, ncol = data$som$grid$xdim)
  plot
}

#' Plot a self-organising map
#'
#' @description Different properties of the units are shown on the map.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param type property to be shown:
#' \itemize{
#'   \item \code{"counts"}: number of samples mapped to each unit
#'   \item \code{"umatrix"}: unified distance matrix, mean Euclidean distance between the prototype of a unit and the prototypes of its direct neighbours. High values separate groups of similar units.
#'   \item \code{"quality"}: mean distance between the samples of a unit and its prototype, see \code{\link{eemsom_samples}}
#'   \item \code{"cluster"}: clusters of the units
#'   \item \code{"samples"}: samples as points on their units, the units are coloured by cluster if available
#' }
#' @param colour optional vector with one value per sample (e.g. sampling site from metadata) used to colour the samples with \code{type = "samples"}
#' @param labels logical, whether sample names are shown with \code{type = "samples"}
#' @param unit_labels logical, whether the numbers of the units are shown. They are used by \code{\link{ggeem.eemsom}}, \code{\link{eemsom_codes}} and \code{\link{eemsom_samples}}.
#'
#' @details Neighbours are defined like in \code{\link[kohonen]{unit.distances}}: six units on hexagonal maps, eight units including the diagonals on rectangular maps.
#'
#' With \code{type = "samples"}, samples of the same unit are arranged in a circle around the centre of the unit.
#'
#' Ejarque-Gonzalez and Butturini (2014) used the U-matrix to find regions of similar units and projected the samples, coloured by sampling site and hydrological conditions, onto the map to explore the patterns among samples.
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @references Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @family plots
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100, seed = 42)
#' eemsom_plot_map(som, unit_labels = TRUE)
#' eemsom_plot_map(som, type = "umatrix")
#' eemsom_plot_map(som, type = "samples", colour = substr(eem_names(eem_list), 1, 2), labels = TRUE)
eemsom_plot_map <- function(som, type = "counts", colour = NULL, labels = FALSE, unit_labels = FALSE){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  type <- match.arg(type, c("counts", "umatrix", "quality", "cluster", "samples"))
  grid <- som$som$grid
  nu <- nrow(grid$pts)
  unit <- som$som$unit.classif
  if(type == "cluster" && is.null(som$clusters)) stop("The units have not been clustered yet!")
  value <- switch(type,
                  counts = tabulate(unit, nbins = nu),
                  umatrix = eemsom_umatrix(som),
                  quality = sapply(seq_len(nu), function(u) mean(sqrt(som$som$distances[unit == u]))),
                  cluster = factor(som$clusters),
                  samples = if(is.null(som$clusters)) rep(NA, nu) else factor(som$clusters))
  poly <- eemsom_polygons(grid)
  poly$value <- value[poly$unit]
  plot <- ggplot(poly, aes(x = x, y = y)) +
    geom_polygon(aes(group = unit, fill = value), colour = "grey50") +
    coord_equal() +
    theme_void()
  if(is.numeric(value)){
    plot <- plot +
      scale_fill_viridis_c(na.value = "white") +
      labs(fill = switch(type, counts = "samples", umatrix = "distance", quality = "distance"))
  } else if(type == "samples" && is.null(som$clusters)){
    plot <- plot +
      scale_fill_manual(values = "white", na.value = "white", guide = "none")
  } else {
    plot <- plot +
      scale_fill_viridis_d() +
      labs(fill = "cluster")
  }
  if(type == "samples"){
    smp <- eemsom_samplepos(som)
    if(!is.null(colour)){
      if(length(colour) != nrow(smp)) stop("colour needs one value per sample!")
      smp$colour <- colour
      plot <- plot + geom_point(data = smp, aes(colour = colour)) + labs(colour = NULL)
      # without the light yellow end of viridis, points are visible on the white background
      plot <- plot + if(is.numeric(colour)) scale_colour_viridis_c(end = 0.9) else scale_colour_viridis_d(end = 0.9)
    } else {
      plot <- plot + geom_point(data = smp)
    }
    if(labels) plot <- plot + geom_text(data = smp, aes(label = sample), size = 2.5, vjust = -0.8)
  }
  if(unit_labels){
    # unit numbers in the upper part of the units, so they do not overlap with the samples
    units <- data.frame(unit = seq_len(nu), x = grid$pts[, 1], y = grid$pts[, 2] + 0.38)
    plot <- plot + geom_text(data = units, aes(label = unit), size = 3, colour = "grey30")
  }
  plot
}

#' Plot component planes of a self-organising map
#'
#' @description A component plane shows the value of one variable in the prototypes of all units. Similar component planes indicate correlated variables. For maps of EEMs, the variables are wavelength pairs, for maps of PARAFAC loadings (\code{\link{eempf_som}}), the variables are the PARAFAC components.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}} or \code{\link{eempf_som}}
#' @param ex numeric vector of excitation wavelengths (maps of EEMs only)
#' @param em numeric vector of emission wavelengths, same length as \code{ex}. Each pair of \code{ex} and \code{em} gives one component plane.
#' @param comps numeric vector or names of PARAFAC components (maps of PARAFAC loadings only). NULL (default) plots all components.
#' @param scale \code{"common"} (default) uses one colour scale for all component planes, so values can be compared between planes. \code{"each"} scales each component plane from 0 (minimum) to 1 (maximum), which shows patterns of weak peaks or components more clearly.
#'
#' @details The closest wavelengths available in the data are used. With \code{scale = "common"}, values are on the scale of the normalised data.
#'
#' Ejarque-Gonzalez and Butturini (2014) compared the component planes of all wavelength pairs by their correlation to discriminate fluorescence components, see \code{\link{eemsom_regions}}.
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @references Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @family plots
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100, seed = 42)
#' # fluorescence peaks B, T, A and C (Coble 1996)
#' eemsom_plot_cplane(som, ex = c(275, 275, 260, 350), em = c(310, 340, 450, 450))
#' eemsom_plot_cplane(som, ex = c(275, 275, 260, 350), em = c(310, 340, 450, 450), scale = "each")
#'
#' data(pf_models)
#' pfsom <- eempf_som(pf4[[1]], rlen = 100, seed = 42)
#' eemsom_plot_cplane(pfsom, scale = "each")
eemsom_plot_cplane <- function(som, ex = NULL, em = NULL, comps = NULL, scale = "common"){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  scale <- match.arg(scale, c("common", "each"))
  codes <- som$som$codes[[1]]
  poly <- eemsom_polygons(som$som$grid)
  if(inherits(som, "eempf_som")){
    if(!is.null(ex) || !is.null(em)) stop("Please use comps instead of ex and em for maps of PARAFAC loadings!")
    if(is.null(comps)) comps <- colnames(codes)
    if(is.numeric(comps)) comps <- colnames(codes)[comps]
    if(any(is.na(comps)) || !all(comps %in% colnames(codes))) stop("comps need to be numbers or names of the PARAFAC components!")
    vars <- comps
    labels <- comps
    fill_lab <- if(som$normalise) "relative\nFmax" else "Fmax"
  } else {
    if(!is.null(comps)) stop("comps can only be used for maps of PARAFAC loadings, please use ex and em!")
    if(is.null(ex) || is.null(em) || length(ex) != length(em)) stop("ex and em need to have the same length!")
    exs <- sapply(ex, function(w) som$ex[which.min(abs(som$ex - w))])
    ems <- sapply(em, function(w) som$em[which.min(abs(som$em - w))])
    vars <- paste0(exs, "|", ems)
    removed <- !vars %in% colnames(codes)
    if(any(removed)) stop("The wavelength pair ex = ", exs[removed][1], " nm, em = ", ems[removed][1], " nm was removed due to missing values!")
    labels <- paste0("ex ", exs, " nm, em ", ems, " nm")
    fill_lab <- "fluorescence\n(normalised)"
  }
  table <- lapply(seq_along(vars), function(i){
    value <- unname(codes[, vars[i]])
    if(scale == "each"){
      # a constant plane is shown in the middle of the colour scale
      value <- if(diff(range(value)) > 0) (value - min(value)) / diff(range(value)) else rep(0.5, length(value))
    }
    poly %>%
      mutate(value = value[unit], pair = labels[i])
  }) %>%
    bind_rows() %>%
    mutate(pair = factor(pair, levels = unique(labels)))
  ggplot(table, aes(x = x, y = y)) +
    geom_polygon(aes(group = unit, fill = value), colour = "grey50") +
    scale_fill_viridis_c() +
    coord_equal() +
    theme_void() +
    facet_wrap(~ pair) +
    labs(fill = if(scale == "each") "relative\n(per plane)" else fill_lab)
}

#' Plot the training progress of a self-organising map
#'
#' @description The mean distance between samples and their closest unit during the training is plotted. The curve should flatten at the end of the training, otherwise \code{rlen} might be increased.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @family plots
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100, seed = 42)
#' eemsom_plot_training(som)
eemsom_plot_training <- function(som){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  data.frame(iteration = seq_len(nrow(som$som$changes)), change = som$som$changes[, 1]) %>%
    ggplot(aes(x = iteration, y = change)) +
    geom_line() +
    theme_minimal() +
    labs(x = "iteration", y = "mean distance to closest unit")
}

#' Polygons of the units of a map
#'
#' @param grid object of class somgrid
#'
#' @return data frame with columns unit, x and y
#' @noRd
eemsom_polygons <- function(grid){
  pts <- grid$pts
  if(grid$topo == "hexagonal"){
    ang <- seq(30, 330, by = 60) * pi / 180
    dx <- cos(ang) / sqrt(3)
    dy <- sin(ang) / sqrt(3)
  } else {
    dx <- c(-0.5, 0.5, 0.5, -0.5)
    dy <- c(-0.5, -0.5, 0.5, 0.5)
  }
  data.frame(unit = rep(seq_len(nrow(pts)), each = length(dx)),
             x = rep(pts[, 1], each = length(dx)) + dx,
             y = rep(pts[, 2], each = length(dx)) + dy)
}

#' Unified distance matrix of a map
#'
#' @param som object of class eemsom
#'
#' @return numeric vector, mean distance of each prototype to the prototypes of the neighbouring units
#' @importFrom stats dist
#' @noRd
eemsom_umatrix <- function(som){
  codes <- som$som$codes[[1]]
  neighbours <- eemsom_unitdist(som$som$grid) == 1
  cdist <- as.matrix(dist(codes))
  sapply(seq_len(nrow(codes)), function(u) mean(cdist[u, neighbours[u, ]]))
}

#' Positions of samples on the map
#'
#' @param som object of class eemsom
#'
#' @return data frame with columns sample, unit, x and y
#' @noRd
eemsom_samplepos <- function(som){
  smp <- eemsom_samples(som)
  for(u in unique(smp$unit)){
    idx <- which(smp$unit == u)
    if(length(idx) > 1){
      ang <- 2 * pi * (seq_along(idx) - 1) / length(idx)
      smp$x[idx] <- smp$x[idx] + 0.25 * cos(ang)
      smp$y[idx] <- smp$y[idx] + 0.25 * sin(ang)
    }
  }
  smp[, c("sample", "unit", "x", "y")]
}
