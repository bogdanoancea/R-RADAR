# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea
# Native R spectral records and file dispatch.

.io_components <- c("baseline", "cosmic_rays", "noise", "peaks")
.io_or <- function(x, y) if (is.null(x)) y else x
.io_need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("This format requires the R package '", pkg, "'. Install it with install.packages('", pkg, "').", call. = FALSE)
}
.io_matrix <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (length(unique(lengths(x))) > 1L) return(lapply(x, as.numeric))
    return(do.call(rbind, lapply(x, as.numeric)))
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (ncol(x) == 1L && is.null(dim(attr(x, "original")))) x else x
}
.io_rows <- function(x) if (is.list(x)) x else lapply(seq_len(nrow(x)), function(i) x[i, ])
.io_array <- function(x) {
  if (is.null(dim(x))) return(matrix(as.numeric(x), nrow = 1L))
  if (length(dim(x)) <= 2L) return(.io_matrix(x))
  d <- dim(x)
  # Last axis is spectral, preceding axes follow R's normal array order.
  matrix(as.numeric(x), ncol = d[length(d)])
}

#' Create a spectral data record
#' @param intensities Numeric vector, spectra-by-channel matrix, array with the
#'   spectral axis last, or a list of spectra of unequal lengths.
#' @param wavenumbers Spectral axis, or a list of axes for unequal lengths.
#' @param names One name per spectrum.
#' @param positions Matrix of spatial or temporal coordinates.
#' @param position_names Names of coordinate columns.
#' @param map_shape Dimensions of the measurement grid.
#' @param x_units,y_units Axis and intensity units.
#' @param meta Named list of metadata.
#' @param arrays Optional named list of multiple spectral arrays.
#' @param source Source file name.
#' @return A spectral_data S3 object; each matrix has one spectrum per row.
#' @export
#' @examples
#' set.seed(1)
#'
#' # Two spectra sharing one axis.
#' obj <- spectral_data(matrix(runif(2 * 16), nrow = 2),
#'                      wavenumbers = seq(400, 1800, length.out = 16),
#'                      names = c("a", "b"),
#'                      x_units = "cm-1", y_units = "counts")
#' obj
#' dim(obj$intensities)
spectral_data <- function(intensities, wavenumbers = NULL, names = NULL,
                          positions = NULL, position_names = character(),
                          map_shape = NULL, x_units = "", y_units = "", meta = list(),
                          arrays = NULL, source = "") {
  if (missing(intensities)) intensities <- arrays[[1L]]
  if (is.null(arrays)) arrays <- list(intensity = intensities)
  if (is.null(base::names(arrays)) || any(!nzchar(base::names(arrays))))
    stop("arrays must be a named list", call. = FALSE)
  d <- dim(intensities)
  if (is.null(map_shape) && length(d) > 2L) map_shape <- d[-length(d)]
  arrays <- lapply(arrays, function(a) if (is.list(a)) .io_matrix(a) else .io_array(a))
  first <- .io_or(arrays$intensity, .io_or(arrays$spectra, arrays[[1L]]))
  n <- if (is.list(first)) length(first) else nrow(first)
  lens <- if (is.list(first)) lengths(first) else rep(ncol(first), n)
  if (!is.null(wavenumbers)) {
    if (is.list(wavenumbers)) {
      if (length(wavenumbers) != n || any(lengths(wavenumbers) != lens))
        stop("wavenumber axes must match the spectrum lengths", call. = FALSE)
    } else if (length(wavenumbers) != lens[1L] || any(lens != length(wavenumbers)))
      stop("wavenumbers must match the spectrum length", call. = FALSE)
  }
  if (!is.null(names) && length(names) != n) stop("names must match the number of spectra", call. = FALSE)
  if (!is.null(positions)) {
    positions <- as.matrix(positions)
    if (nrow(positions) != n && length(positions) == n) positions <- matrix(positions, ncol=1L)
    if (nrow(positions) != n) stop("positions must have one row per spectrum", call. = FALSE)
    if (!length(position_names)) position_names <- .io_or(colnames(positions), paste0("position", seq_len(ncol(positions))))
  }
  if (!is.null(map_shape) && prod(map_shape) != n) stop("prod(map_shape) must equal the number of spectra", call. = FALSE)
  if (is.null(map_shape) && !is.null(positions) && ncol(positions) >= 2L) {
    nx <- length(unique(positions[,1L])); ny <- length(unique(positions[,2L]))
    if (nx > 1L && ny > 1L && nx*ny == n && nrow(unique(positions[,1:2,drop=FALSE])) == n) {
      xx <- matrix(positions[,1L], nrow=ny, byrow=TRUE)
      yy <- matrix(positions[,2L], nrow=ny, byrow=TRUE)
      if (all(xx == rep(xx[1L,], each=ny)) && all(yy == yy[,1L])) map_shape <- c(ny,nx)
    }
  }
  if(!is.null(map_shape)&&is.null(meta$map_order))meta$map_order<-"F"
  structure(list(arrays=arrays, spectra=first, intensities=first,
    wavenumbers=wavenumbers, names=if(is.null(names)) NULL else as.character(names),
    positions=positions, position_names=as.character(position_names), map_shape=map_shape,
    x_units=x_units, y_units=y_units, meta=meta, metadata=meta, source=source),
    class=c("spectral_data", "radar_spectral_data"))
}

#' @rdname spectral_data
#' @export
SpectralData <- spectral_data

#' @export
print.spectral_data <- function(x, ...) {
  a <- x$intensities
  cat("SpectralData:", if(is.list(a)) length(a) else nrow(a), "spectra;",
      if(is.list(a)) "unequal lengths" else paste(ncol(a), "channels"), "\n")
  if (nzchar(x$source)) cat("Source:", x$source, "\n")
  invisible(x)
}

.io_record <- function(obj, wavenumbers = NULL, names = NULL) {
  if (inherits(obj, "spectral_data")) return(obj)
  if (is.list(obj) && !is.null(obj$arrays)) {
    return(spectral_data(arrays=obj$arrays, wavenumbers=.io_or(wavenumbers,obj$wavenumbers),
      names=.io_or(names,obj$names), positions=obj$positions,
      position_names=.io_or(obj$position_names,character()), map_shape=obj$map_shape,
      x_units=.io_or(obj$x_units,""), y_units=.io_or(obj$y_units,""), meta=.io_or(obj$meta,list())))
  }
  if (is.list(obj) && all(.io_components %in% base::names(obj))) {
    arrays <- c(list(spectra=.io_or(obj$spectra,obj$X)), unclass(obj)[.io_components])
    meta <- .io_or(obj$meta,list())
    if (!is.null(obj$params)) meta <- c(meta,list(params=obj$params,kind="synthetic_dataset"))
    if (!is.null(obj$model)) meta$model <- obj$model
    if (!is.null(obj$units)) meta$units <- obj$units
    return(spectral_data(arrays=arrays, wavenumbers=.io_or(wavenumbers,obj$wavenumbers),
      names=.io_or(names,obj$names), positions=obj$positions,
      position_names=.io_or(obj$position_names,character()), map_shape=obj$map_shape,
      x_units=.io_or(obj$x_units,""), y_units=.io_or(obj$units,""), meta=meta))
  }
  spectral_data(obj, wavenumbers=wavenumbers, names=names)
}

.io_formats <- list(
  csv=list(ext=c("csv"), pkg="", description="Delimited text; wide, rows or long layout"),
  tsv=list(ext=c("tsv","tab"), pkg="", description="Tab-separated text"),
  text=list(ext=c("txt","dat","asc","prn","xy"), pkg="", description="Instrument text exports"),
  rruff=list(ext="rruff", pkg="", description="RRUFF text (read only)"),
  rds=list(ext="rds", pkg="", description="Native R spectral record"),
  json=list(ext="json",pkg="jsonlite",description="JSON spectral records"),
  npy=list(ext="npy",pkg="",description="NumPy numeric array, native R reader/writer"),
  npz=list(ext="npz",pkg="jsonlite",description="NumPy archive, native R reader/writer"),
  spc=list(ext="spc",pkg="",description="Galactic SPC, new and old formats"),
  jcamp=list(ext=c("jdx","dx","jcamp","jcm"),pkg="",description="JCAMP-DX AFFN, SQZ, DIF and DUP"),
  hdf5=list(ext=c("h5","hdf5","he5","nc","nc4"),pkg="hdf5r, jsonlite",description="HDF5 / NetCDF-4"),
  sqlite=list(ext=c("db","sqlite","sqlite3"),pkg="DBI, RSQLite, jsonlite",description="SQLite spectra and metadata"),
  excel=list(ext=c("xlsx","xlsm","xls"),pkg="readxl (read); openxlsx, jsonlite (write)",description="Excel workbook; writing requires openxlsx"),
  parquet=list(ext=c("parquet","pq"),pkg="arrow, jsonlite",description="Apache Parquet"),
  feather=list(ext=c("feather","arrow"),pkg="arrow, jsonlite",description="Arrow IPC file (Feather V2); Feather V1 is read"),
  matlab=list(ext="mat",pkg="jsonlite",description="MATLAB level 5 numeric arrays and v7.3 HDF5"),
  wdf=list(ext="wdf",pkg="",description="Renishaw WiRE binary spectra and axes"),
  labspec=list(ext="xml",pkg="xml2",description="Horiba LabSpec XML"),
  trivista=list(ext="tvf",pkg="xml2",description="Princeton TriVista XML")
)

#' List supported spectral file formats
#' @param show Print the format table.
#' @return A data frame with extensions, capabilities and optional R packages.
#' @export
#' @examples
#' # Every adapter, with the optional package each one needs.
#' caps <- formats()
#' head(caps)
#'
#' # Printed summary instead of a returned value.
#' formats(show = TRUE)
formats <- function(show = FALSE) {
  out <- do.call(rbind,lapply(base::names(.io_formats),function(n) {
    f <- .io_formats[[n]]
    data.frame(format=n,extensions=paste0(".",f$ext,collapse=" "),read=TRUE,
      write=!n %in% c("rruff","wdf","labspec","trivista"),requires=if(nzchar(f$pkg)) f$pkg else "-",
      description=f$description, stringsAsFactors=FALSE)
  }))
  rownames(out) <- NULL
  if(show) print(out,row.names=FALSE)
  out
}
.io_format <- function(path, format=NULL) {
  key <- tolower(sub("^\\.","",.io_or(format,tools::file_ext(path))))
  if (key %in% base::names(.io_formats)) return(key)
  for (n in base::names(.io_formats)) if(key %in% .io_formats[[n]]$ext) return(n)
  stop("Unknown spectral format: '",key,"'. See formats().",call.=FALSE)
}

#' Read spectra from a file
#' @param path Input path.
#' @param format Format name or extension; inferred when omitted.
#' @param ... Format options including layout, delimiter, decimal, skip_header,
#'   encoding, sheet, key, table and query.
#' @return A spectral_data record.
#' @export
#' @examples
#' path <- system.file("extdata", "synthetic_spectra.csv", package = "radar.raman")
#' data <- read(path)
#' data
#' dim(data$intensities)
#' head(data$wavenumbers)
read <- function(path, format=NULL, ...) {
  if(!file.exists(path)) stop("File does not exist: ",path,call.=FALSE)
  fmt <- .io_format(path,format)
  out <- switch(fmt,
    csv=.io_read_text(path,...),tsv=.io_read_text(path,default_delimiter="\t",...),
    text=.io_read_text(path,...),rruff=.io_read_text(path,...),
    rds=.io_record(readRDS(path)),json=.io_read_json(path,...),
    npy=.io_interpret(list(data=.io_npy_read(path))),npz=.io_read_npz(path,...),
    spc=.io_read_spc(path,...),jcamp=.io_read_jcamp(path,...),
    hdf5=.io_read_h5(path,...),sqlite=.io_read_sqlite(path,...),excel=.io_read_excel(path,...),
    parquet=.io_read_arrow(path,feather=FALSE,...),feather=.io_read_arrow(path,feather=TRUE,...),
    matlab=.io_read_mat(path,...),wdf=.io_read_wdf(path,...),
    labspec=.io_read_vendor_xml(path,"labspec",...),trivista=.io_read_vendor_xml(path,"trivista",...))
  out$intensities <- out$spectra <- .io_or(out$arrays$intensity,.io_or(out$arrays$spectra,out$arrays[[1L]]))
  out$source <- normalizePath(path,mustWork=FALSE)
  if(is.null(out$meta$format)) out$meta$format <- fmt
  out$metadata <- out$meta
  out
}

#' Write spectra or a decomposition
#' @param obj Spectral data, components, synthetic dataset or numeric array.
#' @param path Output path.
#' @param format Format name or extension.
#' @param wavenumbers,names Optional spectral axis and names for a plain array.
#' @param ... Format options.
#' @return The path, or paths when a one-array format splits the components.
#' @export
#' @examples
#' path <- system.file("extdata", "synthetic_spectra.csv", package = "radar.raman")
#' data <- read(path)
#'
#' # Round trip through a core format; no optional packages required.
#' out <- file.path(tempdir(), "spectra.jdx")
#' write(data, out)
#' back <- read(out)
#' dim(back$intensities)
#' unlink(out)
write <- function(obj, path, format=NULL, wavenumbers=NULL, names=NULL, ...) {
  fmt <- .io_format(path,format)
  if(fmt %in% c("rruff","wdf","labspec","trivista")) stop(fmt," is a read-only format",call.=FALSE)
  rec <- .io_record(obj,wavenumbers,names)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(fmt %in% c("spc","jcamp","npy") && length(rec$arrays)>1L) {
    paths <- vapply(base::names(rec$arrays),function(n) {
      p <- file.path(dirname(path),paste0(tools::file_path_sans_ext(basename(path)),"_",n,".",tools::file_ext(path)))
      sub <- rec; sub$arrays <- setNames(list(rec$arrays[[n]]),"intensity")
      sub$intensities <- sub$spectra <- sub$arrays[[1L]]
      write(sub,p,format=fmt,...)
    }, character(1))
    return(unname(paths))
  }
  switch(fmt,csv=.io_write_text(path,rec,delimiter=",",...),
    tsv=.io_write_text(path,rec,delimiter="\t",...),text=.io_write_text(path,rec,delimiter="\t",...),
    rds=saveRDS(rec,path),json=.io_write_json(path,rec,...),
    npy=.io_npy_write(path,rec$arrays[[1L]]),npz=.io_write_npz(path,rec,...),
    spc=.io_write_spc(path,rec,...),jcamp=.io_write_jcamp(path,rec,...),
    hdf5=.io_write_h5(path,rec,...),sqlite=.io_write_sqlite(path,rec,...),
    excel=.io_write_excel(path,rec,...),parquet=.io_write_arrow(path,rec,FALSE,...),
    feather=.io_write_arrow(path,rec,TRUE,...),matlab=.io_write_mat(path,rec,...))
  path
}

#' @rdname read
#' @export
read_spectra <- read
#' Save spectra to a file
#' @param path Output path.
#' @param spectra Spectra-by-channel matrix.
#' @param wavenumbers,names Optional axis and names.
#' @param ... Format options.
#' @return Output path(s).
#' @export
save_spectra <- function(path,spectra,wavenumbers=NULL,names=NULL,...) write(spectra,path,wavenumbers=wavenumbers,names=names,...)
#' Save components to a file
#' @param components Components object.
#' @param path Output path.
#' @param ... Format options.
#' @return Output path(s).
#' @export
save_components <- function(components,path,...) write(components,path,...)
#' Read saved components
#' @param path Input path; one component path for single-array formats.
#' @param format Optional format override.
#' @param ... Format options.
#' @return A radar_components object.
#' @export
#' @examples
#' dataset <- generate(2, length = 128, seed = 7)
#' parts <- decompose(dataset$X, model = "standard")
#'
#' path <- file.path(tempdir(), "components.rds")
#' save_components(parts, path)
#' restored <- load_components(path)
#' isTRUE(all.equal(restored$peaks, parts$peaks))
#' unlink(path)
load_components <- function(path,format=NULL,...) {
  fmt <- .io_format(path,format)
  if(fmt %in% c("spc","jcamp","npy")) {
    stem <- sub("_(spectra|baseline|cosmic_rays|noise|peaks)$","",tools::file_path_sans_ext(path))
    arr <- lapply(c("spectra",.io_components),function(k) read(paste0(stem,"_",k,".",tools::file_ext(path)),format=fmt,...) )
    rec <- arr[[1L]]; rec$arrays <- setNames(lapply(arr,`[[`,"intensities"),c("spectra",.io_components))
  } else rec <- read(path,format=fmt,...)
  if(!all(c("spectra",.io_components) %in% base::names(rec$arrays))) stop("File does not contain all four components and spectra",call.=FALSE)
  do.call(radar_components,c(rec$arrays[c("spectra",.io_components)],list(wavenumbers=rec$wavenumbers,
    units=.io_or(rec$meta$units,"original"),model=.io_or(rec$meta$model,""),meta=rec$meta,
    names=rec$names,positions=rec$positions,position_names=rec$position_names,map_shape=rec$map_shape,x_units=rec$x_units)))
}
#' Save a synthetic dataset
#' @param dataset Synthetic dataset.
#' @param path Output path; .npz is appended when no extension is present.
#' @param ... Format options.
#' @return Output path(s).
#' @export
save_dataset <- function(dataset,path,...) {if(!nzchar(tools::file_ext(path))) path<-paste0(path,".npz");write(dataset,path,...)}
#' Load a synthetic dataset
#' @param path Input path.
#' @param format Optional format override.
#' @param ... Format options.
#' @return A synthetic_dataset object.
#' @export
load_dataset <- function(path,format=NULL,...) {
  if(!nzchar(tools::file_ext(path))) path<-paste0(path,".npz")
  rec <- read(path,format,...)
  x <- .io_or(rec$arrays$spectra,rec$arrays$intensity)
  if(is.null(x)||!all(.io_components %in% base::names(rec$arrays))) stop("File is not a RADAR synthetic dataset",call.=FALSE)
  do.call(synthetic_dataset,c(list(X=x),rec$arrays[.io_components],list(params=.io_or(rec$meta$params,rec$meta))))
}
