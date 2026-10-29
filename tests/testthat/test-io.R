test_that("io checks pass", {
  local({
    folder <- tempfile("radar-io-"); dir.create(folder)
    on.exit(unlink(folder, recursive=TRUE))
    X <- matrix(sin(seq_len(60)/4), 6, 10)
    wn <- c(100,101,103,106,110,115,121,128,136,145)
    coords <- cbind(x=rep(0:2,2),y=rep(0:1,each=3))
    d <- spectral_data(X,wn,names=paste0("point",1:6),positions=coords,
      position_names=c("x","y"),map_shape=c(2,3),x_units="cm-1",y_units="counts",
      meta=list(instrument="test",note="accent: é"))
    has_json <- requireNamespace("jsonlite",quietly=TRUE)
    extensions <- c("csv","tsv","txt","npy","spc","jdx","rds")
    if(has_json)extensions<-c(extensions,"json","npz","mat")
    optional <- list(h5=c("hdf5r","jsonlite"), db=c("DBI","RSQLite","jsonlite"),
      xlsx=c("readxl","openxlsx","jsonlite"), parquet=c("arrow","jsonlite"), feather=c("arrow","jsonlite"))
    for(ext in names(optional)) if(all(vapply(optional[[ext]],requireNamespace,logical(1),quietly=TRUE))) extensions<-c(extensions,ext)
    for(ext in extensions) {
      path <- file.path(folder,paste0("spectra.",ext))
      write(d,path)
      back <- read(path)
      tol <- if(ext=="spc")1e-6 else 1e-10
      expect_true(all(identical(dim(back$intensities),dim(X)),max(abs(back$intensities-X))<tol))
      if(ext!="npy")expect_true(all(max(abs(back$wavenumbers-wn))<tol))
      if(ext%in%c("json","npz","h5","db","xlsx","mat","rds","parquet","feather"))
        expect_true(all(identical(as.integer(back$map_shape),c(2L,3L)),back$x_units=="cm-1",
          back$meta$instrument=="test",identical(back$names,d$names)))
    }
    # Wide tables with more spectra than channels have an explicit spectral axis.
    short <- spectral_data(matrix(seq_len(30),6,5),c(100,101,105,110,120))
    path <- file.path(folder,"wide.csv"); write(short,path,layout="wide")
    expect_true(all(identical(dim(read(path)$spectra),c(6L,5L)),all(read(path)$spectra==short$spectra)))
    # Plain decimal commas and instrument preambles remain readable.
    path <- file.path(folder,"european.csv")
    writeLines(c("Acquisition export","Raman shift;Intensity","100,5;1,25","101,5;2,25","102,5;3,25"),path)
    back <- read(path)
    expect_true(all(all(back$wavenumbers==c(100.5,101.5,102.5)),all(back$spectra==c(1.25,2.25,3.25))))
    # Rich formats preserve all training targets and generation parameters.
    ds <- generate(3,length=32,seed=4)
    for(ext in c("rds",if(has_json)c("json","npz","mat"))) {
      path <- file.path(folder,paste0("dataset.",ext));save_dataset(ds,path)
      back <- load_dataset(path)
      expect_true(all(max(abs(back$X-ds$X))<1e-12,max(abs(back$noise-ds$noise))<1e-12,
        back$params$length==32))
    }
    # Decompositions preserve all four components; one-array formats split files.
    components <- radar_components(ds$X,ds$baseline,ds$cosmic_rays,ds$noise,ds$peaks,
      wavenumbers=seq_len(32),model="standard",units="original",meta=list(note="known targets"))
    for(ext in extensions) {
      path<-file.path(folder,paste0("components.",ext));paths<-save_components(components,path)
      back<-load_components(paths[1L])
      tol<-if(ext=="spc")1e-6 else 1e-10
      for(k in c("spectra","baseline","cosmic_rays","noise","peaks"))
        expect_true(all(identical(dim(back[[k]]),dim(components[[k]])),max(abs(back[[k]]-components[[k]]))<tol))
      if(ext%in%c("json","npz","h5","db","xlsx","mat","rds","parquet","feather"))
        expect_true(all(back$model=="standard",back$units=="original",back$meta$note=="known targets"))
    }
    # Unequal-length spectra round-trip without resampling or truncation.
    ragged <- spectral_data(list(c(1,2,3),c(4,5,6,7)),list(c(10,20,30),c(40,50,60,70)))
    for(ext in c("csv","jdx","spc","rds",if(has_json)c("json","npz"))) {
      path<-file.path(folder,paste0("ragged.",ext));write(ragged,path);back<-read(path)
      expect_true(all(is.list(back$spectra),identical(lengths(back$spectra),c(3L,4L)),
        all(back$spectra[[1L]]==c(1,2,3)),all(back$spectra[[2L]]==c(4,5,6,7))))
    }
    if(requireNamespace("xml2",quietly=TRUE)) {
      path<-file.path(folder,"labspec.xml")
      writeLines(c('<LSX_Data>','<LSX_Tree><LSX ID="0x7B697861"><LSX ID="0x1">',
        '<LSX ID="0x7D6CD4DB">100 200 300 400</LSX><LSX ID="0x7C696E75">1/cm</LSX>',
        '</LSX></LSX></LSX_Tree><LSX_Matrix>',
        '<LSX_Row Size="4">1 2 3 4</LSX_Row><LSX_Row Size="4">5 6 7 8</LSX_Row>',
        '</LSX_Matrix></LSX_Data>'),path)
      back<-read(path);expect_true(all(all(back$wavenumbers==c(100,200,300,400)),all(back$spectra==rbind(1:4,5:8))))
      path<-file.path(folder,"trivista.tvf")
      writeLines(c('<XmlMain><Documents><Document Label="Intensity" DataLabel="Counts">',
        '<xDim><Axis ValueArray="4|100|200|300|400" Unit="1/cm" /></xDim>',
        '<Data><Frame TimeStamp="10000000">1;2;3;4</Frame><Frame TimeStamp="20000000">5;6;7;8</Frame></Data>',
        '</Document></Documents></XmlMain>'),path)
      back<-read(path);expect_true(all(all(back$wavenumbers==c(100,200,300,400)),all(back$spectra==rbind(1:4,5:8)),all(back$positions[,1L]==c(0,1))))
    }
  })
})

