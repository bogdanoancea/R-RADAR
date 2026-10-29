# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjöberg; Copyright (c) 2026 Bogdan Oancea
# MATLAB level 5 files are parsed directly, including zlib compressed elements.
.io_mat_tag <- function(buf,pos,endian="little") {
  pair<-.io_bin(buf,pos,2L,2L,signed=FALSE,endian=endian)
  if(endian=="big")pair<-rev(pair)
  if(pair[2L]>0L&&pair[2L]<=4L&&pair[1L]>0L&&pair[1L]<=18L)
    return(list(type=pair[1L],data=buf[seq.int(pos+5L,pos+4L+pair[2L])],next_pos=pos+8L))
  kind<-.io_bin(buf,pos,endian=endian);n<-.io_bin(buf,pos+4L,endian=endian)
  if(n<0L||pos+8L+n>length(buf))stop("Invalid MATLAB data element",call.=FALSE)
  dat<-if(n)buf[seq.int(pos+9L,pos+8L+n)]else raw()
  # Compressed elements do not carry outer alignment bytes.
  list(type=kind,data=dat,next_pos=pos+8L+if(kind==15L)n else ceiling(n/8)*8)
}
.io_mat_num <- function(tag,endian="little") {
  sizes<-c(`1`=1L,`2`=1L,`3`=2L,`4`=2L,`5`=4L,`6`=4L,`7`=4L,`9`=8L,`12`=8L,`13`=8L,`16`=1L,`17`=2L,`18`=4L)
  size<-sizes[as.character(tag$type)];if(is.na(size))return(NULL)
  if(!length(tag$data))return(numeric())
  if(tag$type%in%c(5L,6L,12L,13L))return(.io_int_values(tag$data,size,signed=tag$type%in%c(5L,12L),endian=endian))
  .io_bin(tag$data,0L,length(tag$data)/size,size,
    if(tag$type%in%c(7L,9L))"double"else"integer",signed=!tag$type%in%c(2L,4L,6L,16L,17L,18L)||size>=4L,endian=endian)
}
.io_mat_matrix <- function(buf,endian="little") {
  pos<-0L;get<-function(){t<-.io_mat_tag(buf,pos,endian);pos<<-t$next_pos;t}
  flags<-.io_mat_num(get(),endian);class<-bitwAnd(as.integer(flags[1L]),255L)
  if(bitwAnd(as.integer(flags[1L]),2048L))stop("Complex MATLAB spectra are not supported",call.=FALSE)
  dims<-as.integer(.io_mat_num(get(),endian));n<-if(length(dims))prod(dims)else 0L
  nm<-rawToChar(get()$data)
  if(class==2L) {
    size<-as.integer(.io_mat_num(get(),endian))[1L];namesraw<-get()$data
    fn<-vapply(seq_len(length(namesraw)%/%size),function(i).io_text(namesraw,(i-1L)*size,size),character(1))
    values<-vector("list",n*length(fn))
    for(i in seq_along(values)){tag<-get();values[[i]]<-if(tag$type==14L).io_mat_matrix(tag$data,endian)$value else NULL}
    ans<-if(n==1L)setNames(values,fn)else lapply(seq_len(n),function(i)setNames(values[seq.int((i-1L)*length(fn)+1L,i*length(fn))],fn))
    return(list(name=nm,value=ans))
  }
  if(class==1L) {
    vals<-vector("list",n)
    for(i in seq_len(n)){tag<-get();vals[[i]]<-if(tag$type==14L).io_mat_matrix(tag$data,endian)$value else NULL}
    return(list(name=nm,value=vals))
  }
  if(pos>=length(buf))return(list(name=nm,value=numeric()))
  real<-get();value<-.io_mat_num(real,endian)
  if(class==4L) {
    if(real$type==16L)value<-rawToChar(real$data)
    else if(is.null(value))value<-""
    else if(length(dims)>=2L&&dims[1L]>1L){v<-matrix(as.integer(value),nrow=dims[1L]);value<-apply(v,1L,function(v)intToUtf8(v[v!=0L]))}
    else value<-intToUtf8(as.integer(value[value!=0L]))
  }else if(!is.null(value)&&length(dims)&&prod(dims)==length(value))dim(value)<-dims
  list(name=nm,value=value)
}
.io_mat_flatten <- function(x,prefix="") {
  if(is.numeric(x)||is.logical(x))return(setNames(list(x),sub("/$","",prefix)))
  out<-list()
  if(is.list(x))for(k in seq_along(x))out<-c(out,.io_mat_flatten(x[[k]],paste0(prefix,if(is.null(names(x)))k else names(x)[k],"/")))
  out
}
.io_read_mat <- function(path,key=NULL,...) {
  raw<-.io_raw(path)
  if((length(raw)>=8L&&identical(raw[1:8],as.raw(c(137,72,68,70,13,10,26,10)))) ||
     (length(raw)>=520L&&identical(raw[513:520],as.raw(c(137,72,68,70,13,10,26,10))))) return(.io_read_h5(path,key,...))
  if(length(raw)<128L)stop("Not a MATLAB level 5 file",call.=FALSE)
  magic<-rawToChar(raw[127:128]);endian<-if(magic=="IM")"little"else if(magic=="MI")"big"else stop("MATLAB level 4 files are not supported",call.=FALSE)
  values<-list()
  parse<-function(buf,pos=0L) {
    while(pos+8L<=length(buf)) {
      tag<-.io_mat_tag(buf,pos,endian);pos<-tag$next_pos
      if(tag$type==15L)parse(memDecompress(tag$data,"gzip"))
      else if(tag$type==14L){item<-.io_mat_matrix(tag$data,endian);values[[item$name]]<<-item$value}
    }
  }
  parse(raw,128L)
  if(!is.null(values$radar_info)) {
    .io_need("jsonlite");m<-jsonlite::fromJSON(values$radar_info,simplifyVector=FALSE);ns<-as.character(unlist(m$arrays));arrays<-values[ns]
    if(!is.null(key)){if(is.null(arrays[[key]]))stop("Unknown array key",call.=FALSE);arrays<-list(intensity=arrays[[key]])}
    rec<-spectral_data(arrays=arrays,wavenumbers=if(is.null(values$wavenumbers))NULL else as.numeric(values$wavenumbers),names=if(is.null(values$names))NULL else trimws(as.character(unlist(values$names))),positions=values$positions)
    return(.io_apply_meta(rec,m))
  }
  .io_interpret(.io_mat_flatten(values),key)
}
.io_mat_element <- function(type,data) {
  con<-rawConnection(raw(),"wb");writeBin(as.integer(c(type,length(data))),con,size=4L,endian="little");writeBin(data,con)
  if(length(data)%%8L)writeBin(raw(8L-length(data)%%8L),con)
  ans<-rawConnectionValue(con);close(con);ans
}
.io_mat_raw <- function(x,size) {con<-rawConnection(raw(),"wb");writeBin(x,con,size=size,endian="little");v<-rawConnectionValue(con);close(con);v}
.io_mat_encode <- function(name,a) {
  char<-is.character(a);class<-if(char)4L else 6L
  if(char){text<-paste(a,collapse="\n");v<-utf8ToInt(enc2utf8(text));d<-c(1L,length(v));real<-.io_mat_element(4L,.io_mat_raw(as.integer(v),2L))}
  else{d<-dim(a);if(is.null(d))d<-c(1L,length(a));real<-.io_mat_element(9L,.io_mat_raw(as.double(a),8L))}
  body<-c(.io_mat_element(6L,.io_mat_raw(as.integer(c(class,0L)),4L)),.io_mat_element(5L,.io_mat_raw(as.integer(d),4L)),.io_mat_element(1L,charToRaw(name)),real)
  .io_mat_element(14L,body)
}
.io_write_mat <- function(path,rec,...) {
  if(any(vapply(rec$arrays,is.list,logical(1))))stop("MATLAB export requires equal-length spectra; use HDF5 or NPZ",call.=FALSE)
  header<-raw(128L);description<-charToRaw("MATLAB 5.0 MAT-file, Platform: R, Created by RADAR native R")
  header[1:116]<-charToRaw(strrep(" ",116L));header[seq_along(description)]<-description
  header[125:128]<-as.raw(c(0L,1L,73L,77L))
  con<-file(path,"wb");on.exit(close(con));writeBin(header,con)
  arrays<-rec$arrays
  if(!is.null(rec$wavenumbers))arrays$wavenumbers<-rec$wavenumbers
  if(!is.null(rec$positions))arrays$positions<-rec$positions
  arrays$radar_info<-.io_meta_json(rec)
  for(k in names(arrays))writeBin(.io_mat_encode(k,arrays[[k]]),con)
  if(!is.null(rec$names))writeBin(.io_mat_cellstr("names",rec$names),con)
}

.io_mat_cellstr <- function(name,values) {
  body<-c(.io_mat_element(6L,.io_mat_raw(as.integer(c(1L,0L)),4L)),
    .io_mat_element(5L,.io_mat_raw(as.integer(c(1L,length(values))),4L)),
    .io_mat_element(1L,charToRaw(name)))
  for(v in values)body<-c(body,.io_mat_encode("",as.character(v)))
  .io_mat_element(14L,body)
}
