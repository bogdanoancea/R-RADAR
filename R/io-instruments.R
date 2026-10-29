# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea
.io_raw <- function(path) readBin(path,"raw",n=file.info(path)$size)
.io_bin <- function(buf,off,n=1L,size=4L,what="integer",signed=TRUE,endian="little") {
  bytes<-n*size
  if(off<0||off+bytes>length(buf))stop("Truncated or invalid binary spectrum file",call.=FALSE)
  con<-rawConnection(buf[seq.int(off+1L,off+bytes)],"rb");on.exit(close(con))
  readBin(con,what,n=n,size=size,signed=signed,endian=endian)
}
.io_text <- function(buf,off,n) {
  r<-buf[seq.int(off+1L,min(length(buf),off+n))];z<-which(r==as.raw(0));if(length(z))r<-head(r,z[1L]-1L)
  trimws(iconv(rawToChar(r),from="latin1",to="UTF-8",sub=""))
}
.io_put <- function(buf,off,value,size=4L,what=NULL,endian="little") {
  if(is.raw(value))raw<-value else {con<-rawConnection(raw(),"wb");writeBin(value,con,size=size,endian=endian);raw<-rawConnectionValue(con);close(con)}
  if(length(raw))buf[seq.int(off+1L,off+length(raw))]<-raw
  buf
}
.io_common_xy <- function(xs,ys,...) {
  if(!length(ys))stop("File contains no spectra",call.=FALSE)
  same<-all(vapply(xs,function(x)isTRUE(all.equal(x,xs[[1L]],tolerance=1e-7,check.attributes=FALSE)),logical(1)))
  spectral_data(if(same)do.call(rbind,ys)else ys,if(same)xs[[1L]]else xs,...)
}
.io_spc_y <- function(buf,pos,n,exp,flags,endian="little") {
  size<-if(exp!= -128L&&bitwAnd(flags,1L))2L else 4L
  y<-if(exp== -128L).io_bin(buf,pos,n,size,"double",endian=endian)else {
    if(pos+n*size>length(buf))stop("Truncated SPC spectrum",call.=FALSE)
    .io_int_values(buf[seq.int(pos+1L,pos+n*size)],size,signed=TRUE,endian=endian)
  }
  if(exp!= -128L)y<-as.numeric(y)*2^(exp-if(size==2L)16L else 32L)
  list(y=y,pos=pos+n*size)
}
.io_read_spc <- function(path,...) {
  buf<-.io_raw(path);if(length(buf)<2L)stop("Not an SPC file",call.=FALSE)
  flags<-as.integer(buf[1L]);version<-as.integer(buf[2L])
  if(!version%in%c(75L,76L,77L))stop("Unknown SPC version",call.=FALSE)
  if(version==77L) {
    exp<-.io_bin(buf,2L,size=2L);n<-as.integer(round(.io_bin(buf,4L,size=4L,what="double")))
    x<-seq(.io_bin(buf,8L,size=4L,what="double"),.io_bin(buf,12L,size=4L,what="double"),length.out=n)
    pos<-224L;ys<-xs<-list();z<-numeric()
    while(pos+32L+4L*n<=length(buf)) {
      subexp<-.io_bin(buf,pos+1L,size=1L);t<-.io_bin(buf,pos+4L,size=4L,what="double");pos<-pos+32L
      e<-if(subexp>0L&&subexp<128L)subexp else exp
      if(e%in%c(128L,-128L))y<-.io_bin(buf,pos,n,4L,"double")
      else {
        r<-matrix(buf[seq.int(pos+1L,pos+4L*n)],nrow=4L)
        r<-as.raw(r[c(3L,4L,1L,2L),,drop=FALSE]);y<-.io_int_values(r,4L,signed=TRUE)*2^(e-32L)
      }
      pos<-pos+4L*n;ys[[length(ys)+1L]]<-y;xs[[length(xs)+1L]]<-x;z<-c(z,t)
      if(!bitwAnd(flags,4L))break
    }
    return(.io_common_xy(xs,ys,positions=if(length(z)>1L)matrix(z,ncol=1L)else NULL,position_names=if(length(z)>1L)"z"else character(),meta=list(format="SPC (old)",comment=.io_text(buf,64L,130L))))
  }
  endian<-if(version==76L)"big"else"little"
  get<-function(off,n=1L,size=4L,what="integer",signed=TRUE).io_bin(buf,off,n,size,what,signed,endian)
  exp<-get(3L,size=1L);n<-get(4L);first<-get(8L,size=8L,what="double");last<-get(16L,size=8L,what="double")
  nsub<-if(bitwAnd(flags,4L))max(1L,get(24L))else 1L;pos<-512L;common<-NULL
  if(!bitwAnd(flags,64L)) {
    if(bitwAnd(flags,128L)){common<-get(pos,n,4L,"double");pos<-pos+4L*n}else common<-seq(first,last,length.out=n)
  }
  xs<-ys<-list();z<-numeric()
  for(i in seq_len(nsub)) {
    if(pos+32L>length(buf))stop("SPC file is missing a subfile header",call.=FALSE)
    subexp<-get(pos+1L,size=1L);t<-get(pos+4L,size=4L,what="double");subn<-get(pos+16L);pos<-pos+32L
    if(bitwAnd(flags,64L)){npts<-subn;x<-get(pos,npts,4L,"double");pos<-pos+4L*npts}else{npts<-n;x<-common}
    value<-.io_spc_y(buf,pos,npts,if(bitwAnd(flags,4L))subexp else exp,flags,endian)
    pos<-value$pos;xs[[i]]<-x;ys[[i]]<-value$y;z<-c(z,t)
  }
  units<-c("arbitrary","cm-1","um","nm","s","min","Hz","kHz","MHz","m/z","ppm","days","years","Raman shift (cm-1)","eV")
  xtype<-get(28L,size=1L,signed=FALSE)
  meta<-list(format="SPC (new)",comment=.io_text(buf,88L,130L),source_instrument=.io_text(buf,45L,9L))
  logoff<-get(248L)
  if(logoff>0L&&logoff+20L<length(buf)) {
    sz<-get(logoff);tx<-get(logoff+8L);if(tx>=20L&&sz>tx&&logoff+tx<length(buf))meta$log<-.io_text(buf,logoff+tx,min(sz-tx,length(buf)-logoff-tx))
  }
  .io_common_xy(xs,ys,positions=if(nsub>1L)matrix(z,ncol=1L)else NULL,position_names=if(nsub>1L)"z"else character(),
    x_units=if(xtype<length(units))units[xtype+1L]else"",y_units="arbitrary intensity",meta=meta)
}
.io_write_spc <- function(path,rec,x_type=13L,y_type=0L,comment="",...) {
  ys<-.io_rows(rec$arrays[[1L]]);n<-length(ys)
  xs<-if(is.list(rec$wavenumbers))rec$wavenumbers else if(!is.null(rec$wavenumbers))rep(list(rec$wavenumbers),n)else lapply(ys,function(y)seq_along(y)-1L)
  common<-all(vapply(xs,function(x)identical(x,xs[[1L]]),logical(1)));x0<-xs[[1L]]
  even<-common&&length(x0)>1L&&max(abs(diff(x0)-(tail(x0,1L)-x0[1L])/(length(x0)-1L)))<1e-6
  flags<-if(n>1L)4L else 0L
  if(!common)flags<-bitwOr(flags,196L)else if(!even)flags<-bitwOr(flags,128L)
  if(n>1L&&!is.null(rec$positions))flags<-bitwOr(flags,16L)
  h<-raw(512L);h[1:4]<-as.raw(c(flags,75L,0L,128L))
  h<-.io_put(h,4L,as.integer(if(common)length(x0)else 0L));h<-.io_put(h,8L,as.double(x0[1L]),8L);h<-.io_put(h,16L,as.double(tail(x0,1L)),8L)
  h<-.io_put(h,24L,as.integer(n));h[29:30]<-as.raw(c(x_type,y_type));h<-.io_put(h,45L,charToRaw("RADAR"))
  if(nzchar(comment))h<-.io_put(h,88L,head(charToRaw(iconv(comment,to="latin1",sub="?")),129L))
  con<-file(path,"wb");on.exit(close(con));writeBin(h,con)
  if(common&&!even)writeBin(as.double(x0),con,size=4L,endian="little")
  for(i in seq_along(ys)) {
    sub<-raw(32L);sub[2L]<-as.raw(128L);sub<-.io_put(sub,2L,as.integer(i-1L),2L)
    sub<-.io_put(sub,4L,as.double(if(is.null(rec$positions))i-1L else rec$positions[i,1L]),4L)
    sub<-.io_put(sub,16L,as.integer(if(common)0L else length(ys[[i]])))
    writeBin(sub,con);if(!common)writeBin(as.double(xs[[i]]),con,size=4L,endian="little")
    writeBin(as.double(ys[[i]]),con,size=4L,endian="little")
  }
}
.io_jcamp_line <- function(line) {
  line<-sub("\\$\\$.*$","",line)
  num<-"[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?|\\?"
  plain<-gsub(num,"",line,perl=TRUE)
  if(!nzchar(gsub("[[:space:],;]","",plain))) {
    ts<-regmatches(line,gregexpr(num,line,perl=TRUE))[[1L]]
    return(list(values=suppressWarnings(as.numeric(ts)),dif=FALSE))
  }
  pat<-"[@A-Ia-i%J-Rj-rS-Zs][0-9]*(?:\\.[0-9]*)?|[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)|\\?"
  ts<-regmatches(line,gregexpr(pat,line,perl=TRUE))[[1L]]
  sqz<-c(setNames(1:9,LETTERS[1:9]),setNames(-(1:9),letters[1:9]),`@`=0L)
  dif<-c(setNames(1:9,LETTERS[10:18]),setNames(-(1:9),letters[10:18]),`%`=0L)
  dup<-c(setNames(1:8,LETTERS[19:26]),s=9L)
  vals<-numeric();kind<-"abs";delta<-0
  decode<-function(v,t)as.numeric(paste0(if(v<0)"-"else"",abs(v),substring(t,2L)))
  for(t in ts) {
    ch<-substring(t,1L,1L)
    if(ch%in%names(sqz)){vals<-c(vals,decode(sqz[ch],t));kind<-"abs"}
    else if(ch%in%names(dif)){delta<-decode(dif[ch],t);vals<-c(vals,tail(vals,1L)+delta);kind<-"dif"}
    else if(ch%in%names(dup)) {
      count<-as.integer(paste0(dup[ch],substring(t,2L)))-1L
      if(count>0L)vals<-c(vals,tail(vals,1L)+if(kind=="dif")delta*seq_len(count)else rep(0,count))
    }else{vals<-c(vals,suppressWarnings(as.numeric(t)));kind<-"abs"}
  }
  list(values=vals,dif=kind=="dif")
}
.io_read_jcamp <- function(path,...) {
  text<-readLines(path,warn=FALSE,encoding="latin1")
  if(any(grepl("^##\\s*NTUPLES",text,ignore.case=TRUE)))stop("JCAMP NTUPLES is not supported",call.=FALSE)
  blocks<-list();stack<-list()
  for(line in text) {
    if(grepl("^\\s*##[^=]*=",line)) {
      label<-toupper(gsub("[[:space:]_/-]","",sub("^\\s*##([^=]*)=.*$","\\1",line)))
      value<-trimws(sub("^[^=]*=","",line))
      if(label=="TITLE") {
        parent<-if(length(stack))stack[[length(stack)]]$labels else list()
        stack[[length(stack)+1L]]<-list(labels=c(parent,list(TITLE=value)),kind=NULL,lines=character());next
      }
      if(!length(stack))stack<-list(list(labels=list(),kind=NULL,lines=character()))
      j<-length(stack)
      if(label=="END") {
        b<-stack[[j]];if(!is.null(b$kind))blocks[[length(blocks)+1L]]<-b
        stack<-head(stack,-1L);next
      }
      stack[[j]]$labels[[label]]<-value
      if(label%in%c("XYDATA","XYPOINTS","PEAKTABLE","DATATABLE")){stack[[j]]$kind<-label;stack[[j]]$form<-value;stack[[j]]$lines<-character()}
    } else if(length(stack)&&!is.null(stack[[length(stack)]]$kind))stack[[length(stack)]]$lines<-c(stack[[length(stack)]]$lines,line)
  }
  for(b in stack)if(!is.null(b$kind))blocks[[length(blocks)+1L]]<-b
  if(!length(blocks))stop("JCAMP file has no spectral data blocks",call.=FALSE)
  xs<-ys<-list();nms<-character()
  for(b in blocks) {
    lab<-b$labels
    number<-function(k,def=NULL){v<-suppressWarnings(as.numeric(strsplit(.io_or(lab[[k]],""),"[[:space:]]+")[[1L]][1L]));if(length(v)&&is.finite(v))v else def}
    xf<-number("XFACTOR",1);yf<-number("YFACTOR",1)
    pairs<-grepl("XY\\.\\.XY|XYXY",gsub(" ","",toupper(b$form)))
    if((b$kind=="XYDATA"&&!pairs)||(b$kind=="DATATABLE"&&grepl("X++",b$form,fixed=TRUE))) {
      y<-numeric();prev<-FALSE;firstline<-NULL
      for(line in b$lines){v<-.io_jcamp_line(line);if(length(v$values)<2L)next;if(is.null(firstline))firstline<-v$values;yy<-v$values[-1L];if(prev&&length(y))yy<-yy[-1L];y<-c(y,yy);prev<-v$dif}
      count<-number("NPOINTS",length(y));if(length(y)>count)y<-head(y,count)
      first<-number("FIRSTX");last<-number("LASTX")
      if(is.null(first)||is.null(last)){first<-firstline[1L]*xf;last<-first+(length(y)-1L)*number("DELTAX",xf)}
      x<-seq(first,last,length.out=length(y));y<-y*yf
    }else {
      num<-"[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?|\\?"
      vals<-unlist(lapply(b$lines,function(l)suppressWarnings(as.numeric(regmatches(sub("\\$\\$.*$","",l),gregexpr(num,sub("\\$\\$.*$","",l),perl=TRUE))[[1L]]))))
      vals<-head(vals,length(vals)%/%2L*2L);xy<-matrix(vals,ncol=2L,byrow=TRUE);x<-xy[,1L]*xf;y<-xy[,2L]*yf
    }
    xs[[length(xs)+1L]]<-x;ys[[length(ys)+1L]]<-y;nms<-c(nms,.io_or(lab$TITLE,paste0("spectrum_",length(xs))))
  }
  lab<-blocks[[1L]]$labels
  .io_common_xy(xs,ys,names=nms,x_units=.io_or(lab$XUNITS,""),y_units=.io_or(lab$YUNITS,""),meta=list(format="JCAMP-DX",labels=lab,jcamp_version=.io_or(lab$JCAMPDX,"")))
}
.io_jfmt <- function(x) ifelse(is.finite(x),sprintf("%.17g",x),"?")
.io_write_jcamp <- function(path,rec,data_type="RAMAN SPECTRUM",...) {
  ys<-.io_rows(rec$arrays[[1L]]);n<-length(ys)
  xs<-if(is.list(rec$wavenumbers))rec$wavenumbers else if(!is.null(rec$wavenumbers))rep(list(rec$wavenumbers),n)else lapply(ys,function(y)seq_along(y)-1L)
  nms<-.io_or(rec$names,paste0("spectrum_",seq_len(n)))
  lines<-if(n>1L)c(paste0("##TITLE=",basename(path)),"##JCAMP-DX=5.01","##DATA TYPE=LINK",paste0("##BLOCKS=",n))else character()
  for(i in seq_len(n)) {
    x<-xs[[i]];y<-ys[[i]];len<-length(y)
    block<-c(paste0("##TITLE=",nms[i]),"##JCAMP-DX=5.01",paste0("##DATA TYPE=",data_type),"##ORIGIN=RADAR",paste0("##XUNITS=",if(nzchar(rec$x_units))toupper(rec$x_units)else"1/CM"),paste0("##YUNITS=",if(nzchar(rec$y_units))rec$y_units else"ARBITRARY UNITS"),"##XFACTOR=1","##YFACTOR=1",paste0("##FIRSTX=",.io_jfmt(x[1L])),paste0("##LASTX=",.io_jfmt(tail(x,1L))),paste0("##NPOINTS=",len))
    even<-len>1L&&max(abs(diff(x)-(tail(x,1L)-x[1L])/(len-1L)))<1e-9
    if(even){block<-c(block,"##XYDATA=(X++(Y..Y))");for(j in seq.int(1L,len,by=8L))block<-c(block,paste(.io_jfmt(c(x[j],y[seq.int(j,min(j+7L,len))])),collapse=" "))}
    else block<-c(block,"##XYPOINTS=(XY..XY)",paste(.io_jfmt(x),.io_jfmt(y),sep=", "))
    lines<-c(lines,block,"##END=")
  }
  if(n>1L)lines<-c(lines,"##END=")
  writeLines(lines,path,useBytes=TRUE)
}
