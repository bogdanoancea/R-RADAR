# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjöberg; Copyright (c) 2026 Bogdan Oancea
.io_meta <- function(rec) list(radar_format=2L,x_units=rec$x_units,y_units=rec$y_units,
  position_names=I(.io_or(rec$position_names,character())),map_shape=if(is.null(rec$map_shape))NULL else I(rec$map_shape),arrays=I(names(rec$arrays)),map_order=.io_or(rec$meta$map_order,"F"),meta=rec$meta,names=if(is.null(rec$names))NULL else I(rec$names))
.io_meta_json <- function(rec) { .io_need("jsonlite"); as.character(jsonlite::toJSON(.io_meta(rec),auto_unbox=TRUE,null="null",digits=NA)) }
.io_apply_meta <- function(rec,js) {
  .io_need("jsonlite")
  m<-if(is.list(js))js else tryCatch(jsonlite::fromJSON(js,simplifyVector=FALSE),error=function(e)list())
  if(!is.list(m))return(rec)
  for(k in c("x_units","y_units"))if(!is.null(m[[k]]))rec[[k]]<-as.character(m[[k]])
  if(!is.null(m$position_names))rec$position_names<-as.character(unlist(m$position_names,use.names=FALSE))
  if(!is.null(m$map_shape))rec$map_shape<-as.integer(unlist(m$map_shape))
  if(!is.null(m$names)&&is.null(rec$names))rec$names<-as.character(unlist(m$names))
  if(!is.null(m$meta))rec$meta<-utils::modifyList(.io_or(rec$meta,list()),m$meta)
  if(!is.null(rec$map_shape))rec$meta$map_order<-.io_or(m$map_order,.io_or(m$meta$map_order,"C"))
  rec$metadata<-rec$meta;rec
}
.io_numeric_flat <- function(v) {
  if(is.null(v))return(NA_real_)
  if(is.list(v))return(unlist(lapply(v,.io_numeric_flat),use.names=FALSE))
  as.numeric(v)
}
.io_json_array <- function(x) {
  if(is.null(x))return(NULL)
  if(!is.list(x))return(x)
  if(!length(x))return(numeric())
  if(all(vapply(x,function(v)is.atomic(v)&&is.null(dim(v))&&length(v)<=1L,logical(1))))return(.io_numeric_flat(x))
  vals<-lapply(x,.io_numeric_flat)
  if(length(unique(lengths(vals)))==1L)do.call(rbind,vals)else vals
}
.io_interpret <- function(arrays,key=NULL) {
  arrays<-arrays[vapply(arrays,function(a)is.numeric(a)||is.logical(a),logical(1))]
  if(!length(arrays))stop("No numeric spectral arrays found",call.=FALSE)
  leaves<-sub(".*/","",names(arrays));names(leaves)<-names(arrays)
  axis_alias<-c("wavenumbers","wavenumber","raman_shift","ramanshift","shift","axis","axisscale","x","xdata","wavelength","wavelengths","spectral_axis")
  ak<-which(tolower(leaves)%in%axis_alias & vapply(arrays,function(a)is.null(dim(a))||sum(dim(a)>1L)<=1L,logical(1)))
  axis<-if(length(ak))as.numeric(arrays[[ak[1L]]])else NULL
  cmp<-c("spectra",.io_components); have<-cmp[cmp%in%leaves]
  if(!is.null(key)) {
    where<-match(key,names(arrays));if(is.na(where))where<-match(key,leaves)
    if(is.na(where))stop("No array named '",key,"'; available: ",paste(names(arrays),collapse=", "),call.=FALSE)
    out<-list(intensity=arrays[[where]])
  } else if(length(have)>=2L)out<-setNames(lapply(have,function(c)arrays[[match(c,leaves)]]),have)
  else {
    cand<-match(c("intensity","intensities","spectra","spectrum","x_data","X","data","y","counts","signal"),leaves)
    cand<-cand[!is.na(cand)&!cand%in%ak]
    if(!length(cand)){other<-setdiff(seq_along(arrays),ak);if(!length(other))other<-seq_along(arrays);cand<-other[which.max(lengths(arrays[other]))]}
    out<-list(intensity=arrays[[cand[1L]]])
  }
  map<-NULL
  out<-lapply(out,function(a){
    d<-dim(a)
    if(!is.null(axis)&&length(d)>=2L&&tail(d,1L)!=length(axis)&&length(axis)%in%d) {
      channel<-which(d==length(axis))[1L]
      a<-aperm(a,c(setdiff(seq_along(d),channel),channel));d<-dim(a)
    }
    if(is.null(axis)&&length(d)==2L&&nrow(a)>ncol(a)*4L)a<-t(a)
    if(length(dim(a))>2L) {
      d<-dim(a);map<<-head(d,-1L)
      return(matrix(as.vector(aperm(a,c(rev(seq_len(length(d)-1L)),length(d)))),ncol=tail(d,1L)))
    }
    .io_array(a)
  })
  if(!is.null(axis)&&length(axis)!=ncol(out[[1L]]))axis<-NULL
  spectral_data(arrays=out,wavenumbers=axis,map_shape=map,meta=if(is.null(map))list()else list(map_order="C"))
}
.io_read_json <- function(path,key=NULL,...) {
  .io_need("jsonlite")
  obj<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!is.null(obj$radar_format)) {
    arrays<-lapply(obj$arrays,.io_json_array)
    if(!is.null(key)){if(is.null(arrays[[key]]))stop("Unknown array key",call.=FALSE);arrays<-list(intensity=arrays[[key]])}
    wn<-.io_json_array(obj$wavenumbers);if(is.matrix(wn)&&nrow(wn)==1L)wn<-as.numeric(wn)
    rec<-spectral_data(arrays=arrays,wavenumbers=wn,names=if(is.null(obj$names))NULL else as.character(unlist(obj$names)),positions=.io_json_array(obj$positions),
      position_names=as.character(unlist(obj$position_names)),map_shape=if(is.null(obj$map_shape))NULL else as.integer(unlist(obj$map_shape)),
      x_units=.io_or(obj$x_units,""),y_units=.io_or(obj$y_units,""),meta=.io_or(obj$meta,list()))
    return(.io_apply_meta(rec,obj))
  }
  if(is.list(obj)&&is.null(names(obj))&&length(obj)&&is.list(obj[[1L]])&&!is.null(names(obj[[1L]]))) {
    ys<-lapply(obj,function(it){k<-names(it)[tolower(names(it))%in%c("y","intensity","spectra","spectrum","counts","signal")][1L];as.numeric(unlist(it[[k]]))})
    xs<-lapply(seq_along(obj),function(i){it<-obj[[i]];k<-names(it)[tolower(names(it))%in%c("x","wavenumbers","wavenumber","axis","wavelength")];if(length(k))as.numeric(unlist(it[[k[1L]]]))else seq_along(ys[[i]])-1L})
    same<-all(vapply(xs,function(x)isTRUE(all.equal(x,xs[[1L]])),logical(1)))
    return(spectral_data(if(same)do.call(rbind,ys)else ys,if(same)xs[[1L]]else xs,names=vapply(seq_along(obj),function(i)as.character(.io_or(obj[[i]]$name,.io_or(obj[[i]]$title,paste0("spectrum_",i)))),character(1))))
  }
  if(is.list(obj)&&!is.null(names(obj))) .io_interpret(lapply(obj,.io_json_array),key) else .io_interpret(list(data=.io_json_array(obj)),key)
}
.io_write_json <- function(path,rec,...) {
  .io_need("jsonlite")
  obj<-.io_meta(rec);obj$arrays<-lapply(rec$arrays,function(a)if(is.matrix(a))unname(a)else a)
  obj$wavenumbers<-if(is.null(rec$wavenumbers)||is.list(rec$wavenumbers))rec$wavenumbers else I(rec$wavenumbers);obj$positions<-rec$positions
  jsonlite::write_json(obj,path,auto_unbox=TRUE,null="null",digits=NA,na="null",matrix="rowmajor")
}
.io_pack <- function(a) {
  if(!is.list(a))return(list(data=a,lengths=NULL))
  lens<-lengths(a);out<-matrix(NaN,length(a),max(lens))
  for(i in seq_along(a))out[i,seq_len(lens[i])]<-a[[i]]
  list(data=out,lengths=lens)
}
.io_unpack <- function(a,lens=NULL) if(is.null(lens))a else lapply(seq_along(lens),function(i)as.numeric(a[i,seq_len(lens[i])]))

# NumPy format is decoded directly; object/pickle arrays are never evaluated.
.io_npy_read <- function(path) {
  con<-if(inherits(path,"connection"))path else file(path,"rb")
  own<-!inherits(path,"connection");if(own)on.exit(close(con))
  if(!identical(readBin(con,"raw",6L),c(as.raw(147),charToRaw("NUMPY"))))stop("Invalid NumPy magic",call.=FALSE)
  ver<-readBin(con,"integer",2L,size=1L,signed=FALSE)
  size<-if(ver[1L]==1L)2L else 4L
  len<-readBin(con,"integer",1L,size=size,signed=if(size==2L)FALSE else TRUE,endian="little")
  header<-rawToChar(readBin(con,"raw",len))
  dtype<-sub(".*['\"]descr['\"]\\s*:\\s*['\"]([^'\"]+)['\"].*","\\1",header,perl=TRUE)
  shape<-sub(".*['\"]shape['\"]\\s*:\\s*\\(([^)]*)\\).*","\\1",header,perl=TRUE)
  dims<-suppressWarnings(as.integer(trimws(strsplit(shape,",",fixed=TRUE)[[1L]])));dims<-dims[!is.na(dims)]
  n<-if(length(dims))prod(dims)else 1L
  code<-substring(dtype,2L,2L); bytes<-as.integer(substring(dtype,3L));endi<-if(startsWith(dtype,">"))"big"else"little"
  if(code=="O"||is.na(bytes))stop("Object/pickle NumPy arrays are not supported",call.=FALSE)
  if(code=="f")a<-readBin(con,"double",n,size=bytes,endian=endi)
  else if(code%in%c("i","u","b","?")) {
    if(bytes>=4L){r<-readBin(con,"raw",n*bytes);a<-.io_int_values(r,bytes,signed=code=="i",endian=endi)}
    else a<-readBin(con,"integer",n,size=bytes,signed=code=="i"||bytes>=4L,endian=endi)
  } else if(code=="S")a<-vapply(seq_len(n),function(i){r<-readBin(con,"raw",bytes);r<-r[seq_len(if(any(r==as.raw(0)))which(r==as.raw(0))[1L]-1L else length(r))];rawToChar(r)},character(1))
  else if(code=="U")a<-vapply(seq_len(n),function(i){v<-readBin(con,"integer",bytes,size=4L,endian=endi);intToUtf8(v[v!=0L])},character(1))
  else stop("Unsupported NumPy dtype: ",dtype,call.=FALSE)
  if(length(a)!=n)stop("Truncated NumPy array",call.=FALSE)
  if(length(dims)>1L) {
    if(grepl("fortran_order['\"]\\s*:\\s*True",header,perl=TRUE))dim(a)<-dims
    else a<-aperm(array(a,rev(dims)),rev(seq_along(dims)))
  }
  a
}
.io_npy_write <- function(path,a,scalar=FALSE) {
  if(is.list(a))stop("NumPy .npy requires a rectangular numeric array",call.=FALSE)
  con<-file(path,"wb");on.exit(close(con))
  dims<-if(scalar)integer()else dim(a);if(is.null(dims))dims<-length(a)
  shape<-if(length(dims)==1L)paste0(dims,",")else paste(dims,collapse=", ")
  dtype<-if(is.character(a))paste0("<U",max(1L,vapply(a,function(x)length(utf8ToInt(x)),integer(1))))else"<f8"
  h<-paste0("{'descr': '",dtype,"', 'fortran_order': False, 'shape': (",shape,"), }")
  pad<-(16L-((10L+nchar(h,type="bytes")+1L)%%16L))%%16L
  h<-paste0(h,strrep(" ",pad),"\n")
  writeBin(c(as.raw(147),charToRaw("NUMPY"),as.raw(c(1,0))),con)
  writeBin(as.integer(nchar(h,type="bytes")),con,size=2L,endian="little")
  writeBin(charToRaw(h),con)
  vals<-if(length(dims)>1L)as.vector(aperm(a,rev(seq_along(dims))))else as.vector(a)
  if(is.character(a)) {
    width<-as.integer(substring(dtype,3L))
    for(x in vals){v<-utf8ToInt(x);writeBin(as.integer(c(v,rep(0L,width-length(v)))),con,size=4L,endian="little")}
  } else writeBin(as.double(vals),con,size=8L,endian="little")
  invisible(path)
}
.io_read_npz <- function(path,key=NULL,...) {
  listing<-utils::unzip(path,list=TRUE)$Name
  listing<-listing[grepl("\\.npy$",listing)]
  files<-lapply(listing,function(n){con<-unz(path,n,open="rb");on.exit(close(con));.io_npy_read(con)})
  names(files)<-sub("\\.npy$","",listing)
  if(!is.null(files$`_radar`)) {
    .io_need("jsonlite");m<-jsonlite::fromJSON(as.character(files$`_radar`),simplifyVector=FALSE)
    ns<-as.character(unlist(m$arrays));arrays<-setNames(lapply(ns,function(k).io_unpack(files[[paste0("array_",k)]],files[[paste0("lengths_",k)]])),ns)
    if(!is.null(key)){if(is.null(arrays[[key]]))stop("Unknown array key",call.=FALSE);arrays<-list(intensity=arrays[[key]])}
    wn<-.io_unpack(files$wavenumbers,files$wavenumber_lengths)
    rec<-spectral_data(arrays=arrays,wavenumbers=wn,names=files$names,positions=files$positions)
    return(.io_apply_meta(rec,m))
  }
  if(!is.null(files$X)&&all(.io_components%in%names(files))){files$spectra<-files$X;files$X<-NULL}
  rec<-.io_interpret(files,key)
  meta<-.io_or(files$`_meta`,files$`_params`)
  if(!is.null(meta)){.io_need("jsonlite");rec$meta<-tryCatch(jsonlite::fromJSON(as.character(meta)),error=function(e)list())}
  rec
}
.io_write_npz <- function(path,rec,...) {
  folder<-tempfile("radar-npz-");dir.create(folder);on.exit(unlink(folder,recursive=TRUE))
  files<-list(`_radar`=.io_meta_json(rec))
  for(k in names(rec$arrays)){p<-.io_pack(rec$arrays[[k]]);files[[paste0("array_",k)]]<-p$data;if(!is.null(p$lengths))files[[paste0("lengths_",k)]]<-p$lengths}
  if(!is.null(rec$wavenumbers)){p<-.io_pack(rec$wavenumbers);files$wavenumbers<-p$data;if(!is.null(p$lengths))files$wavenumber_lengths<-p$lengths}
  if(!is.null(rec$names))files$names<-rec$names
  if(!is.null(rec$positions))files$positions<-rec$positions
  for(k in names(files)).io_npy_write(file.path(folder,paste0(k,".npy")),files[[k]],scalar=k=="_radar")
  dest<-file.path(normalizePath(dirname(path),mustWork=TRUE),basename(path));old<-setwd(folder);on.exit(setwd(old),add=TRUE)
  if(file.exists(dest))unlink(dest)
  status<-utils::zip(dest,paste0(names(files),".npy"),flags="-q")
  if(status!=0L)stop("Could not write NumPy archive; ensure the system zip utility is available",call.=FALSE)
}
.io_h5_value <- function(ds) {
  a<-ds$read();if(length(dim(a))>1L)aperm(a,rev(seq_along(dim(a))))else a
}
.io_h5_string_type <- function() {
  t<-hdf5r::h5types$H5T_C_S1$copy();t$set_size(Inf);t$set_cset("UTF-8");t
}
.io_h5_set <- function(f,k,a) {
  if(length(dim(a))>1L)a<-aperm(a,rev(seq_along(dim(a))))
  if(is.character(a))f$create_dataset(k,robj=enc2utf8(a),dtype=.io_h5_string_type())else f[[k]]<-a
}
.io_read_h5 <- function(path,key=NULL,...) {
  .io_need("hdf5r");f<-hdf5r::H5File$new(path,"r");on.exit(f$close_all())
  get<-function(k)if(f$exists(k)).io_h5_value(f[[k]])else NULL
  has<-hdf5r::h5attr_names(f)
  if("_radar"%in%has||"radar_format"%in%has) {
    js<-if("_radar"%in%has)hdf5r::h5attr(f,"_radar")else NULL
    .io_need("jsonlite");m<-if(is.null(js))list(arrays=intersect(c("spectra",.io_components,"intensity"),f$ls()$name))else jsonlite::fromJSON(js,simplifyVector=FALSE)
    ns<-as.character(unlist(m$arrays));arrays<-setNames(lapply(ns,function(k).io_unpack(get(k),get(paste0(k,"_lengths")))),ns)
    if(!is.null(key)){if(is.null(arrays[[key]]))stop("Unknown array key",call.=FALSE);arrays<-list(intensity=arrays[[key]])}
    rec<-spectral_data(arrays=arrays,wavenumbers=.io_unpack(get("wavenumbers"),get("wavenumber_lengths")),names=get("names"),positions=get("positions"))
    return(.io_apply_meta(rec,m))
  }
  if(is.null(key)&&f$exists("Regions")) {
    ls<-f[["Regions"]]$ls();parts<-lapply(ls$name,function(n).io_foreign_array(get(paste0("Regions/",n,"/Dataset"))))
    return(spectral_data(do.call(rbind,parts)))
  }
  ls<-f$ls(recursive=TRUE);ns<-ls$name[ls$obj_type=="H5I_DATASET"]
  arrays<-setNames(lapply(ns,get),ns);.io_interpret(arrays,key)
}
.io_write_h5 <- function(path,rec,...) {
  .io_need("hdf5r");f<-hdf5r::H5File$new(path,"w");on.exit(f$close_all())
  for(k in names(rec$arrays)){p<-.io_pack(rec$arrays[[k]]);.io_h5_set(f,k,p$data);if(!is.null(p$lengths))f[[paste0(k,"_lengths")]]<-p$lengths}
  if(!is.null(rec$wavenumbers)){p<-.io_pack(rec$wavenumbers);.io_h5_set(f,"wavenumbers",p$data);if(!is.null(p$lengths))f[["wavenumber_lengths"]]<-p$lengths}
  if(!is.null(rec$names)).io_h5_set(f,"names",rec$names)
  if(!is.null(rec$positions)).io_h5_set(f,"positions",rec$positions)
  f$create_attr("radar_format",robj=2L,space=hdf5r::H5S$new("scalar"))
  f$create_attr("_radar",robj=enc2utf8(.io_meta_json(rec)),dtype=.io_h5_string_type(),space=hdf5r::H5S$new("scalar"))
}

.io_int_values <- function(buf,size,signed=TRUE,endian="little") {
  n<-length(buf)/size
  vapply(seq_len(n),function(i) {
    b<-as.integer(buf[seq.int((i-1L)*size+1L,i*size)])
    if(endian=="big")b<-rev(b)
    if(signed&&b[size]>=128L) -sum((255L-b)*256^(seq_len(size)-1L))-1 else sum(b*256^(seq_len(size)-1L))
  },numeric(1))
}

.io_foreign_array <- function(a) {
  d<-dim(a)
  if(length(d)<=2L)return(.io_array(a))
  matrix(as.vector(aperm(a,c(rev(seq_len(length(d)-1L)),length(d)))),ncol=tail(d,1L))
}
