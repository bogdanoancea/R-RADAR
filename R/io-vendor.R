# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Bogdan Oancea
# Independent native readers of the vendor data grammars. Field identifiers
# are format constants, not executable third-party parsers.
# WDF structure references: Alex Henderson, DOI:10.5281/zenodo.495477;
# https://github.com/alchem0x2A/py-wdf-reader (block/offset documentation).
.io_u64 <- function(buf,off) sum(as.integer(buf[seq.int(off+1L,off+8L)])*256^(0:7))
.io_read_wdf <- function(path,...) {
  buf<-.io_raw(path)
  if(length(buf)<512L||rawToChar(buf[1:4])!="WDF1")stop("Not a Renishaw WDF file",call.=FALSE)
  blocks<-list();pos<-0
  while(pos+16<=length(buf)) {
    tag<-rawToChar(buf[seq.int(pos+1L,pos+4L)]);size<-.io_u64(buf,pos+8L)
    if(!is.finite(size)||size<16||pos+size>length(buf))stop("Invalid WDF block length",call.=FALSE)
    blocks[[tag]]<-c(offset=pos,size=size);pos<-pos+size
  }
  if(is.null(blocks$DATA)||is.null(blocks$XLST))stop("WDF is missing DATA or XLST blocks",call.=FALSE)
  points<-.io_bin(buf,60L);capacity<-.io_u64(buf,64L);count<-.io_u64(buf,72L)
  if(points<1L||count<1||count*points>(blocks$DATA["size"]-16L)/4L)stop("Invalid WDF spectrum dimensions",call.=FALSE)
  spectra<-matrix(.io_bin(buf,blocks$DATA["offset"]+16L,count*points,4L,"double"),nrow=count,byrow=TRUE)
  xoffset<-blocks$XLST["offset"];axis<-.io_bin(buf,xoffset+24L,points,4L,"double")
  unit<-.io_bin(buf,xoffset+20L)
  units<-c("","cm-1","nm","nm","eV","um","counts","electrons","mm","m","K","Pa","s","ms","h","days","pixels","intensity")
  positions<-NULL;pn<-character();origins<-list()
  if(!is.null(blocks$ORGN)) {
    p<-blocks$ORGN["offset"]+20L;num<-.io_bin(buf,blocks$ORGN["offset"]+16L)
    for(i in seq_len(num)) {
      typ<-bitwAnd(.io_bin(buf,p),2147483647L);label<-.io_text(buf,p+8L,16L)
      if(!is.na(typ)&&typ%in%c(3L,4L,5L,11L)) {
        if(typ==11L){a<-vapply(seq_len(count),function(j).io_u64(buf,p+24L+(j-1L)*8L),numeric(1));a<-(a-a[1L])/1e7}
        else a<-.io_bin(buf,p+24L,count,8L,"double")
        nm<-if(typ==11L)"time"else c("x","y","z")[typ-2L]
        origins[[nm]]<-list(label=label,values=a)
        if(count>1L&&diff(range(a))>0){positions<-cbind(positions,a);pn<-c(pn,nm)}
      }
      p<-p+24L+8L*capacity
    }
  }
  map<-NULL
  if(!is.null(blocks$WMAP)) {
    wh<-.io_bin(buf,blocks$WMAP["offset"]+48L,2L)
    if(all(wh>1L)&&prod(wh)==count)map<-rev(wh)
  }
  spectral_data(spectra,axis,positions=positions,position_names=pn,map_shape=map,
    x_units=if(unit>=0L&&unit<length(units))units[unit+1L]else"",
    y_units="intensity",meta=list(map_order="C",format="Renishaw WDF",title=.io_text(buf,240L,272L),
      application=.io_text(buf,96L,24L),measurement_type=.io_bin(buf,132L),
      laser_wavenumber=.io_bin(buf,156L,size=4L,what="double"),origin_labels=lapply(origins,`[[`,"label"),
      incomplete=count<capacity))
}
.io_xml_numbers <- function(s,sep="[[:space:],;|]+") {
  if(is.na(s)||!nzchar(trimws(s)))return(numeric())
  suppressWarnings(as.numeric(strsplit(trimws(s),sep,perl=TRUE)[[1L]]))
}
.io_xml_attr <- function(node,name,default="") {
  if(inherits(node,"xml_missing"))return(default)
  v<-xml2::xml_attr(node,name);if(is.na(v))default else v
}
.io_xml_id <- function(node,id) {
  if(inherits(node,"xml_missing"))return(node)
  xml2::xml_find_first(node,paste0(".//*[@ID='",id,"']"))
}
.io_read_vendor_xml <- function(path,kind,...) {
  .io_need("xml2")
  doc<-xml2::read_xml(path,options=c("NOBLANKS","NONET"))
  if(kind=="labspec") {
    rows<-xml2::xml_find_all(doc,"//LSX_Matrix[1]/LSX_Row")
    if(!length(rows))stop("LabSpec XML has no LSX_Matrix spectral rows",call.=FALSE)
    ys<-lapply(xml2::xml_text(rows),.io_xml_numbers)
    axes<-.io_xml_id(doc,"0x7B697861")
    signal<-xml2::xml_find_first(axes,"./*[@ID='0x1']")
    wn<-.io_xml_numbers(xml2::xml_text(.io_xml_id(signal,"0x7D6CD4DB")))
    if(!length(wn)||any(!is.finite(wn))||length(wn)!=length(ys[[1L]]))stop("LabSpec XML is missing a compatible calibrated spectral axis",call.=FALSE)
    xunits<-xml2::xml_text(.io_xml_id(signal,"0x7C696E75"));if(is.na(xunits))xunits<-""
    nav<-lapply(c("0x2","0x3"),function(id){node<-xml2::xml_find_first(axes,paste0("./*[@ID='",id,"']"));list(values=.io_xml_numbers(xml2::xml_text(.io_xml_id(node,"0x7D6CD4DB"))),name=xml2::xml_text(.io_xml_id(node,"0x6D707974")))})
    nav<-Filter(function(x)length(x$values)>1L&&all(is.finite(x$values)),nav)
    positions<-NULL;pn<-character();map<-NULL
    if(length(nav)==1L&&length(nav[[1L]]$values)==length(ys)){positions<-matrix(nav[[1L]]$values,ncol=1L);pn<-tolower(nav[[1L]]$name)}
    if(length(nav)==2L&&prod(vapply(nav,function(x)length(x$values),integer(1)))==length(ys)) {
      positions<-as.matrix(expand.grid(nav[[2L]]$values,nav[[1L]]$values))[,2:1,drop=FALSE]
      pn<-vapply(nav,function(x)tolower(x$name),character(1));map<-vapply(rev(nav),function(x)length(x$values),integer(1))
    }
    if(any(is.na(pn)|!nzchar(pn)))pn<-head(c("x","y"),length(pn))
    return(spectral_data(do.call(rbind,ys),wn,positions=positions,position_names=pn,map_shape=map,x_units=xunits,
      meta=list(map_order="C",format="Horiba LabSpec XML",title=xml2::xml_text(.io_xml_id(doc,"0x6C7469D9")))))
  }
  docs<-xml2::xml_find_all(doc,"//Documents/Document")
  if(!length(docs))docs<-xml2::xml_find_all(doc,"//Document")
  if(!length(docs))stop("TriVista XML has no Document",call.=FALSE)
  node<-docs[[1L]]
  frames<-xml2::xml_find_all(node,"./Data/*")
  if(!length(frames))stop("TriVista Document has no Data frames; split glued datasets are not supported",call.=FALSE)
  ys<-lapply(xml2::xml_text(frames),.io_xml_numbers)
  ax<-xml2::xml_find_first(node,"./xDim/*[@ValueArray]")
  vals<-.io_xml_numbers(.io_xml_attr(ax,"ValueArray"))
  if(length(vals)>1L&&vals[1L]==length(vals)-1L)vals<-vals[-1L]
  if(length(vals)!=length(ys[[1L]])||any(!is.finite(vals)))stop("TriVista spectral axis does not match Data frame length",call.=FALSE)
  stamps<-suppressWarnings(as.numeric(xml2::xml_attr(frames,"TimeStamp")))
  pos<-NULL;pn<-character();map<-NULL
  if(length(stamps)>1L&&all(is.finite(stamps))){pos<-matrix((stamps-stamps[1L])/1e7,ncol=1L);pn<-"time"}
  serialized<-.io_xml_attr(node,"InfoSerialized")
  if(nzchar(serialized)) {
    info<-tryCatch(xml2::read_xml(serialized,options=c("NOBLANKS","NONET")),error=function(e)NULL)
    if(!is.null(info)) {
      coords<-lapply(c("X-Axis","Y-Axis"),function(nm) {
        g<-xml2::xml_find_first(info,paste0("//*[@Name='",nm,"']"))
        val<-function(k){v<-xml2::xml_find_first(g,paste0(".//*[@Name='",k,"']"));suppressWarnings(as.numeric(.io_xml_attr(v,"Value")))}
        points<-val("Points");from<-val("From");step<-val("Step")
        if(length(points)&&is.finite(points)&&points>0&&is.finite(from)&&is.finite(step))from+seq.int(0,points-1L)*step else numeric()
      })
      if(length(coords[[1L]])>1L&&length(coords[[2L]])>1L&&prod(lengths(coords))==length(ys)) {
        pos<-as.matrix(expand.grid(x=coords[[1L]],y=coords[[2L]]));pn<-c("x","y");map<-rev(lengths(coords))
      }
    }
  }
  spectral_data(do.call(rbind,ys),vals,positions=pos,position_names=pn,map_shape=map,
    x_units=.io_xml_attr(ax,"Unit"),y_units=.io_xml_attr(node,"DataLabel"),
    meta=list(map_order="C",format="TriVista",label=.io_xml_attr(node,"Label"),record_time=.io_xml_attr(node,"RecordTime")))
}
