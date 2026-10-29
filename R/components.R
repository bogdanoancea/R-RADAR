# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea

#' Names and order of RADAR components
#' @export
COMPONENTS <- c("baseline", "cosmic_rays", "noise", "peaks")

#' Construct a Raman spectrum decomposition
#' @param spectra Input spectra, a vector or a matrix with one spectrum per row.
#' @param baseline,cosmic_rays,noise,peaks Estimated additive components with the same shape.
#' @param wavenumbers Optional spectral axis.
#' @param units Intensity units, original or normalized.
#' @param model Name of the model that produced the result.
#' @param meta Additional metadata.
#' @param names Optional spectrum names.
#' @param positions Optional matrix of scan positions.
#' @param position_names Names of scan coordinate columns.
#' @param map_shape Optional dimensions of a regular spatial map.
#' @param x_units Units of the spectral axis.
#' @return A radar_components object. Derived fields include corrected, denoised,
#'   despiked, baseline_corrected, reconstruction and residual.
#' @export
radar_components <- function(spectra, baseline, cosmic_rays, noise, peaks,
                             wavenumbers = NULL, units = "original", model = "", meta = list(),
                             names = NULL, positions = NULL, position_names = character(),
                             map_shape = NULL, x_units = "") {
  s <- .radar_matrix(spectra)
  a <- list(spectra=if(is.null(dim(spectra))) as.numeric(s) else s, baseline=baseline, cosmic_rays=cosmic_rays, noise=noise, peaks=peaks)
  for (c in COMPONENTS) {
    v <- .radar_matrix(a[[c]], c)
    if (!identical(dim(v),dim(s))) stop("all components must have the shape of spectra", call. = FALSE)
    if (is.null(dim(spectra))) a[[c]] <- as.numeric(v) else a[[c]] <- v
  }
  if (!is.null(wavenumbers) && (!is.numeric(wavenumbers) || length(wavenumbers)!=ncol(s) || any(!is.finite(wavenumbers))))
    stop("wavenumbers must match the spectrum length", call. = FALSE)
  if (!is.null(names) && length(names)!=nrow(s)) stop("names must match the number of spectra",call.=FALSE)
  if (!is.null(positions)) {
    positions <- as.matrix(positions)
    if (nrow(positions)!=nrow(s)) stop("positions must match the number of spectra",call.=FALSE)
    if (length(position_names) && length(position_names)!=ncol(positions))
      stop("position_names must match the coordinate columns",call.=FALSE)
  }
  if (!is.null(map_shape) && (any(map_shape<1) || any(map_shape!=as.integer(map_shape)) || prod(map_shape)!=nrow(s)))
    stop("map_shape must contain positive grid dimensions whose product is the number of spectra",call.=FALSE)
  if (!is.null(map_shape) && is.null(meta$map_order))
    meta$map_order <- if(!is.null(meta$source_metadata$map_order)) meta$source_metadata$map_order else "F"
  structure(c(a,list(wavenumbers=wavenumbers,units=units,model=model,meta=meta,names=names,
                     positions=positions,position_names=position_names,map_shape=map_shape,x_units=x_units)),
            class="radar_components")
}

#' @export
`$.radar_components` <- function(x, name) {
  v <- unclass(x)
  switch(name,
    n_spectra=if(is.null(dim(v$spectra))) 1L else nrow(v$spectra),
    length=if(is.null(dim(v$spectra))) length(v$spectra) else ncol(v$spectra),
    corrected=v$peaks,
    denoised=v$spectra-v$noise,
    despiked=v$spectra-v$cosmic_rays,
    baseline_corrected=v$spectra-v$baseline,
    reconstruction=v$baseline+v$cosmic_rays+v$noise+v$peaks,
    residual=v$spectra-v$baseline-v$cosmic_rays-v$noise-v$peaks,
    metadata=v$meta,
    original=v$spectra,
    v[[name]])
}
#' @export
length.radar_components <- function(x) x$n_spectra
#' @export
print.radar_components <- function(x, ...) {
  cat(sprintf("RADAR components: %d spectrum/spectra x %d samples (%s units; %s model)\n",x$n_spectra,x$length,x$units,x$model))
  if(!is.null(x$map_shape)) cat("Map dimensions:",paste(x$map_shape,collapse=" x "),"\n")
  invisible(x)
}
#' @export
`[.radar_components` <- function(x, i, ...) {
  if(missing(i)) return(x)
  vals <- lapply(c("spectra",COMPONENTS),function(c) .radar_matrix(x[[c]])[i,,drop=FALSE])
  names(vals) <- c("spectra",COMPONENTS)
  ns <- nrow(vals$spectra)
  if(!ns) stop("selection contains no spectra",call.=FALSE)
  meta <- x$meta
  for(k in c("offset","scale")) if(length(meta[[k]])==x$n_spectra) meta[[k]] <- meta[[k]][i]
  do.call(radar_components,c(vals,list(wavenumbers=x$wavenumbers,units=x$units,model=x$model,meta=meta,
    names=if(is.null(x$names)) NULL else x$names[i],
    positions=if(is.null(x$positions)) NULL else x$positions[i,,drop=FALSE],
    position_names=x$position_names,x_units=x$x_units)))
}

.radar_map_order <- function(x) {
  n <- prod(x$map_shape)
  marker <- x$meta$map_order
  if (is.null(marker)) marker <- x$meta$source_metadata$map_order
  unchanged <- if (identical(marker,"C")) {
    shape <- as.integer(x$map_shape)
    as.vector(aperm(array(seq_len(n),dim=rev(shape)),rev(seq_along(shape))))
  } else seq_len(n)
  p <- x$positions
  if (length(x$map_shape)!=2L || is.null(p) || !is.matrix(p) || ncol(p)<2L || nrow(p)!=n)
    return(unchanged)
  labels <- tolower(x$position_names)
  ix <- if("x" %in% labels) match("x",labels) else 1L
  iy <- if("y" %in% labels) match("y",labels) else 2L
  xv <- p[,ix]; yv <- p[,iy]
  if (anyNA(xv) || anyNA(yv) || nrow(unique(cbind(xv,yv)))!=n) return(unchanged)
  nx <- length(unique(xv)); ny <- length(unique(yv))
  if (ny==x$map_shape[1L] && nx==x$map_shape[2L]) return(order(xv,yv))
  if (nx==x$map_shape[1L] && ny==x$map_shape[2L]) return(order(yv,xv))
  unchanged
}

#' Reshape spectral components into a spatial map
#' @param x A decomposition or spectral-data object.
#' @param component Name of a component or derived field.
#' @return An array with dimensions c(map_shape, spectrum_length).
#' @export
as_map <- function(x, component = "peaks") {
  if(is.null(x$map_shape)) stop("these spectra do not come from a regular map",call.=FALSE)
  value <- if(inherits(x,"radar_components")) `$.radar_components`(x,component) else x[[component]]
  if(is.null(value)) stop("unknown component: ",component,call.=FALSE)
  values <- .radar_matrix(value)
  array(values[.radar_map_order(x),,drop=FALSE],dim=c(x$map_shape,ncol(values)))
}

#' Integrate a spectral band
#' @param x A decomposition or spectral-data object with wavenumbers.
#' @param lo,hi Lower and upper spectral band limits.
#' @param component Component to integrate.
#' @return Integrated values, or an array of map_shape for a regular map.
#' @export
band_image <- function(x, lo, hi, component = "peaks") {
  if(is.null(x$wavenumbers)) stop("wavenumbers are needed to select a band",call.=FALSE)
  wn <- x$wavenumbers
  sel <- which(wn >= min(lo,hi) & wn <= max(lo,hi))
  if(!length(sel)) stop("no points between band limits",call.=FALSE)
  value <- if(inherits(x,"radar_components")) `$.radar_components`(x,component) else x[[component]]
  if(is.null(value)) stop("unknown component: ",component,call.=FALSE)
  sel <- sel[order(wn[sel])]; xx <- wn[sel]; yy <- .radar_matrix(value)[,sel,drop=FALSE]
  area <- if(length(sel)==1L) rep(0,nrow(yy)) else
    rowSums(sweep((yy[,-1L,drop=FALSE]+yy[,-ncol(yy),drop=FALSE])/2,2L,diff(xx),"*"))
  if(is.null(x$map_shape)) area else array(area[.radar_map_order(x)],dim=x$map_shape)
}

#' Plot a decomposition
#' @param components A radar_components object.
#' @param index One-based index of the spectrum to plot.
#' @param wavenumbers Optional spectral axis overriding the stored axis.
#' @param title Optional plot title.
#' @param ... Additional arguments passed to the graphics plot function.
#' @return Invisibly, the plotted object.
#' @export
plot_components <- function(components,index=1L,wavenumbers=NULL,title=NULL,...) {
  if(index<1L || index>components$n_spectra) stop("index is outside the spectra",call.=FALSE)
  axis <- if(!is.null(wavenumbers)) wavenumbers else components$wavenumbers
  if(is.null(axis)) axis <- seq_len(components$length)
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  graphics::par(mfrow=c(5L,1L),mar=c(2.5,4,1,1),oma=c(1,0,2,0))
  val <- function(c) .radar_matrix(components[[c]])[index,]
  graphics::plot(axis,val("spectra"),type="l",xlab="",ylab="Input",...)
  graphics::lines(axis,val("baseline"),col="#d98a00")
  cols <- c("#d98a00","#d62728","#7f7f7f","#1f77b4")
  for(j in seq_along(COMPONENTS)) {
    graphics::plot(axis,val(COMPONENTS[j]),type="l",col=cols[j],xlab="",ylab=gsub("_"," ",COMPONENTS[j]),...)
  }
  graphics::mtext(if(is.null(components$wavenumbers)&&is.null(wavenumbers)) "Index" else "Raman shift (cm^-1)",side=1,outer=TRUE)
  graphics::mtext(if(is.null(title)) paste0("RADAR decomposition (",components$model," model)") else title,side=3,outer=TRUE)
  invisible(components)
}
#' @export
plot.radar_components <- function(x,...) plot_components(x,...)

#' Compare predicted and true spectral components
#' @param dataset A synthetic_dataset with ground truth.
#' @param components A decomposition of the same spectra.
#' @param index One-based spectrum index.
#' @param ... Additional graphics parameters.
#' @return Invisibly, the predicted components.
#' @export
plot_comparison <- function(dataset,components,index=1L,...) {
  if(index<1L || index>nrow(dataset$X) || index>components$n_spectra) stop("index is outside the spectra",call.=FALSE)
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  graphics::par(mfrow=c(5L,1L),mar=c(2.5,4,1,1))
  graphics::plot(dataset$X[index,],type="l",xlab="",ylab="Input",...)
  for(c in COMPONENTS) {
    t <- dataset[[c]][index,]; p <- .radar_matrix(components[[c]])[index,]
    graphics::plot(t,type="l",xlab="",ylab=gsub("_"," ",c),ylim=range(t,p),...)
    graphics::lines(p,col="#d62728",lty=2)
  }
  invisible(components)
}
