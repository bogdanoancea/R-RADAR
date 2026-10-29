test_that("io interoperability checks pass", {
  fixture <- function(name) system.file("extdata","interoperability",name,package="RADAR",mustWork=TRUE)
  x <- matrix(seq(0,29)/3-2,nrow=6L,byrow=TRUE)
  axis <- c(100,101,105,110,120)
  
  # Independent NumPy buffers exercise row/column order, endian and integer types.
  expect_true(all(max(abs(RADAR::read(fixture("npy-bef8-F.npy"))$intensities-x))<1e-12,
            identical(as.numeric(RADAR::read(fixture("npy-lei8-C.npy"))$intensities),
                      as.numeric(matrix(-10:19,nrow=6L,byrow=TRUE))),
            identical(as.numeric(RADAR::read(fixture("npy-leu4-C.npy"))$intensities),
                      c(0,1,2^31,2^32-1))))
  
  # Native readers decode compressed SciPy matrices and big-endian small tags.
  expect_true(all(max(abs(RADAR::read(fixture("scipy-True.mat"))$intensities-x))<1e-12,
            max(abs(RADAR::read(fixture("big-endian-small.mat"),key="data")$intensities-x))<1e-12,
            identical(as.numeric(RADAR::read(fixture("scipy-uint32.mat"))$intensities),c(0,1,2^31,2^32-1))))
  
  # SPC integer-coded spectra can contain the signed32 minimum (R's NA sentinel).
  for(endian in c("le","be")) {
    sp <- RADAR::read(fixture(paste0("int32-min-",endian,".spc")))
    expect_true(all(identical(as.numeric(sp$intensities),c(-.5,-2^-32,0,.5-2^-32))))
  }
  
  # A compressed JCAMP Y-check must be dropped before appending the second line.
  jcamp <- RADAR::read(fixture("difdup.jdx"))
  expect_true(all(identical(as.numeric(jcamp$intensities),c(5,5.5,6,6.5,6.5)),
            identical(jcamp$wavenumbers,as.numeric(100:104))))
  
  formats <- c("spc","jdx")
  if(requireNamespace("jsonlite",quietly=TRUE)) formats <- c(formats,"json","npz","mat")
  if(requireNamespace("jsonlite",quietly=TRUE)&&requireNamespace("hdf5r",quietly=TRUE)) formats <- c(formats,"h5")
  if(requireNamespace("jsonlite",quietly=TRUE)&&requireNamespace("DBI",quietly=TRUE)&&requireNamespace("RSQLite",quietly=TRUE)) formats <- c(formats,"db")
  for(ext in formats) {
    obj <- RADAR::read(fixture(paste0("fixture-source.",ext)))
    tol <- if(ext=="spc")1e-6 else 1e-12
    expect_true(all(identical(dim(obj$intensities),dim(x)),max(abs(obj$intensities-x))<tol,
              max(abs(obj$wavenumbers-axis))<tol))
    if(ext %in% c("h5","npz","json","db","mat")) {
      expect_true(all(identical(as.integer(obj$map_shape),c(2L,3L)),
                identical(as.character(obj$names),paste0("p",0:5)),
                max(abs(as_map(obj,"spectra")[2,1,]-x[4,]))<tol))
    }
  }
  
  # Source metadata without physical positions declares a C-order map implicitly.
  for(ext in intersect(formats,c("h5","npz","json","mat","db"))) {
    obj <- RADAR::read(fixture(paste0("fixture-map-no-positions.",ext)))
    expect_true(all(is.null(obj$positions),max(abs(obj$intensities-x))<1e-12,
              max(abs(as_map(obj,"spectra")[2,1,]-x[4,]))<1e-12,
              max(abs(as_map(obj,"spectra")[1,2,]-x[2,]))<1e-12))
    if(ext != "mat") {
      tmp <- tempfile(fileext=paste0(".",ext)); RADAR::write(obj,tmp)
      same <- RADAR::read(tmp)
      expect_true(all(max(abs(as_map(same,"spectra")-as_map(obj,"spectra")))<1e-12))
      unlink(tmp)
    }
  }
  
  # Text instruments with different table orientations retain physical positions.
  wire <- RADAR::read(fixture("wire_map.txt"))
  expect_true(all(identical(dim(wire$intensities),c(6L,20L)),
            identical(as.integer(wire$map_shape),c(2L,3L)),
            as_map(wire,"spectra")[2,1,1]==10,
            as_map(wire,"spectra")[1,2,1]==100))
  lab <- RADAR::read(fixture("labspec_map.txt"))
  expect_true(all(identical(dim(lab$intensities),c(4L,10L)),
            identical(as.numeric(lab$wavenumbers),seq(100,1000,length.out=10))))
  
  # Explicit axis headers disambiguate spectra shorter than the batch size.
  short_wide <- tempfile(fileext=".csv")
  RADAR::write(x,short_wide,wavenumbers=axis,names=paste0("p",0:5),layout="wide")
  expect_true(all(max(abs(RADAR::read(short_wide)$intensities-x))<1e-12))
  unlink(short_wide)
  
  # A minimal independently constructed WDF file retains calibrated map axes.
  wdf <- RADAR::read(fixture("native-test.wdf"))
  expect_true(all(identical(dim(wdf$spectra), c(6L, 5L)),
            all(wdf$spectra == matrix(seq_len(30), 6L, 5L, byrow = TRUE)),
            all(wdf$wavenumbers == axis),
            identical(as.integer(wdf$map_shape), c(2L, 3L)),
            all(wdf$positions[, 1L] == rep(0:2, 2)),
            all(wdf$positions[, 2L] == rep(0:1, each = 3))))
})

