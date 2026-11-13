# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea
.io_key <- function(x) tolower(gsub("[[:space:]]+","",x))
.io_monotonic <- function(x) length(x)>=3L && all(is.finite(x)) && (all(diff(x)>0)||all(diff(x)<0))
.io_wn_alias <- c("wavenumber","wavenumbers","ramanshift","raman_shift","shift","wave","#wave","wl","wavelength","x-axis","xaxis","frequency","cm-1","cm^-1","1/cm","energy","pixel","channel","index")
.io_int_alias <- c("intensity","#intensity","int","counts","count","signal","y","y-axis","yaxis","absorbance","spectrum","value","cps")
.io_pos_alias <- c("x","y","z","#x","#y","#z","t","time","#time","position","pos_x","pos_y","pos_z","xpos","ypos","zpos","depth")
.io_id_alias <- c("spectrum","spectrum_id","spectrumid","id","spec","sample","spectrum_index")
.io_name_alias <- c("name","spectrum_name","label","title")
.io_num <- function(x) suppressWarnings(as.numeric(trimws(gsub('^"|"$','',x))))
.io_numeric_row <- function(x) {
  x <- trimws(x); x <- x[nzchar(x)]
  if(!length(x)) return(FALSE)
  sum(!is.na(.io_num(x))|tolower(x)%in%c("nan","na","?"))>=ceiling(length(x)/2)
}
.io_grid <- function(rows) {
  good <- vapply(rows,.io_numeric_row,logical(1))
  runs <- rle(good); ends <- cumsum(runs$lengths)
  len <- runs$lengths; len[!runs$values]<-0L
  if(!any(len)) stop("No numeric spectral data found",call.=FALSE)
  k <- which.max(len); start <- ends[k]-len[k]+1L; end <- ends[k]
  body <- rows[start:end]; width <- max(lengths(body))
  header <- NULL
  if(start>1L && length(rows[[start-1L]])>=width-1L && !.io_numeric_row(rows[[start-1L]])) header <- rows[[start-1L]]
  nms <- if(is.null(header)) paste0("column_",seq_len(width)) else head(c(header,rep("",width)),width)
  nms[!nzchar(nms)] <- paste0("column_",which(!nzchar(nms)))
  cols <- lapply(seq_len(width),function(j) {
    x <- vapply(body,function(r) if(length(r)>=j) r[j] else "",character(1))
    y <- .io_num(x)
    if(all(!is.na(y)|!nzchar(trimws(x))|tolower(x)%in%c("nan","na","?"))) y else trimws(x)
  }); names(cols)<-make.unique(nms)
  list(cols=as.data.frame(cols,check.names=FALSE,stringsAsFactors=FALSE),header=header)
}
.io_period <- function(x) {
  if(length(x)<6L) return(NULL)
  starts <- which(abs(x-x[1L])<1e-8)
  if(length(starts)<2L) return(NULL)
  len <- starts[2L]-1L
  if(len<3L || length(x)%%len || !.io_monotonic(x[seq_len(len)])) return(NULL)
  if(all(abs(x-rep(x[seq_len(len)],length(x)/len))<1e-8)) len else NULL
}
.io_from_groups <- function(cols, wn, values, id=NULL, nm=NULL, pos=character()) {
  axis <- cols[[wn]]; n <- length(axis)
  if(!is.null(id)) keys <- as.character(cols[[id]])
  else if(length(pos)) keys <- apply(cols[pos],1,paste,collapse="\r")
  else {
    per <- .io_period(axis)
    if(!is.null(per)) keys <- rep(seq_len(n/per),each=per)
    else {
      dif <- sign(diff(axis)); direction <- sign(stats::median(dif)); if(!is.finite(direction)||!direction) direction<-1
      keys <- cumsum(c(1L,as.integer(dif!=direction)))
    }
  }
  ids <- unique(keys); idx <- lapply(ids,function(k) which(keys==k))
  axes <- lapply(idx,function(i) axis[i]); same <- all(vapply(axes,function(a) isTRUE(all.equal(a,axes[[1L]],check.attributes=FALSE)),logical(1)))
  arrays <- lapply(values,function(v) {
    rows <- lapply(idx,function(i) as.numeric(cols[[v]][i]))
    if(same) do.call(rbind,rows) else rows
  }); names(arrays)<-.io_key(values)
  if(length(values)==1L && !names(arrays)%in%c("spectra",.io_components)) names(arrays)<-"intensity"
  nms <- if(!is.null(nm)) vapply(idx,function(i) as.character(cols[[nm]][i[1L]]),character(1)) else NULL
  if(is.null(nms)&&!is.null(id)&&!is.numeric(cols[[id]])) nms<-as.character(ids)
  positions <- if(length(pos)) as.matrix(cols[vapply(idx,`[`,integer(1),1L),pos,drop=FALSE]) else NULL
  spectral_data(arrays=arrays,wavenumbers=if(same) axes[[1L]] else axes,names=nms,
    positions=positions,position_names=sub("^#","",.io_key(pos)),meta=list(layout="long"))
}
.io_table <- function(cols,header=names(cols),layout="auto") {
  cols <- as.data.frame(cols,check.names=FALSE,stringsAsFactors=FALSE)
  nums <- names(cols)[vapply(cols,is.numeric,logical(1))]
  if(!length(nums)) stop("Table contains no numeric columns",call.=FALSE)
  key <- .io_key(names(cols))
  find <- function(aliases) {i<-which(key%in%aliases);if(length(i))names(cols)[i[1L]]else NULL}
  wn <- if(!is.null(header)) find(.io_wn_alias) else NULL
  if(layout%in%c("auto","long") && !is.null(wn) && wn%in%nums) {
    vals <- nums[nums!=wn & .io_key(nums)%in%c("spectra",.io_components,setdiff(.io_int_alias,c(.io_id_alias,.io_pos_alias)))]
    if(!length(vals)) vals <- nums[nums!=wn & .io_key(nums)%in%setdiff(.io_int_alias,.io_id_alias)]
    id <- find(.io_id_alias); if(!is.null(id)&&id%in%vals) id<-NULL
    nm <- find(.io_name_alias)
    pos <- nums[.io_key(nums)%in%.io_pos_alias & !nums%in%c(wn,vals)]
    if(length(vals)&&(anyDuplicated(cols[[wn]])||!is.null(id)||layout=="long")) return(.io_from_groups(cols,wn,vals,id,nm,pos))
  }
  if(layout%in%c("auto","long")&&is.null(header)&&length(nums)>=3L&&length(nums)==ncol(cols)) {
    wn <- nums[length(nums)-1L]
    if(!is.null(.io_period(cols[[wn]]))||layout=="long") {
      ans <- .io_from_groups(cols,wn,tail(nums,1L),pos=head(nums,-2L))
      ans$position_names <- if(ncol(cols)==4L)c("y","x") else head(c("x","y","z","t"),ncol(cols)-2L)
      return(ans)
    }
  }
  tab <- as.matrix(cols[nums]); nr<-nrow(tab);nc<-ncol(tab)
  columns <- function() {
    if(nc>=2L&&.io_monotonic(tab[,1L])) spectral_data(t(tab[,-1L,drop=FALSE]),tab[,1L],if(is.null(header))NULL else nums[-1L],meta=list(layout="wide"))
    else spectral_data(t(tab),names=if(is.null(header))NULL else nums,meta=list(layout="wide"))
  }
  rows <- function() {
    lead <- if(all(is.na(tab[1L,])))nc else which(!is.na(tab[1L,]))[1L]-1L
    if(nr>=2L&&lead<nc&&.io_monotonic(tab[1L,seq.int(lead+1L,nc)])) {
      spectral_data(tab[-1L,seq.int(lead+1L,nc),drop=FALSE],tab[1L,seq.int(lead+1L,nc)],
        positions=if(lead)tab[-1L,seq_len(lead),drop=FALSE]else NULL,
        position_names=head(c("x","y","z","t"),lead),meta=list(layout="rows"))
    } else if(!is.null(header)&&all(is.finite(.io_num(nums)))) spectral_data(tab,.io_num(nums),meta=list(layout="rows"))
    else spectral_data(tab,meta=list(layout="rows"))
  }
  if(layout%in%c("wide","columns")) return(columns())
  if(layout=="rows") return(rows())
  if(!layout%in%c("auto","long")) stop("layout must be auto, wide, columns, rows or long",call.=FALSE)
  if(nr>=2L&&is.na(tab[1L,1L]))return(rows())
  if(!is.null(header)&&nc>2L&&.io_monotonic(.io_num(nums)))return(rows())
  if(!is.null(header)&&.io_key(nums[1L])%in%.io_wn_alias&&nc>=2L&&.io_monotonic(tab[,1L]))return(columns())
  if(nc<=2L||nr>nc)columns()else rows()
}
.io_read_text <- function(path,layout="auto",delimiter=NULL,decimal=NULL,comments="#",skip_header=0L,encoding=NULL,default_delimiter=NULL,...) {
  raw <- readBin(path,"raw",n=file.info(path)$size)
  if(length(raw)>=2L&&identical(raw[1:2],as.raw(c(255,254)))) enc <- "UTF-16LE"
  else if(length(raw)>=2L&&identical(raw[1:2],as.raw(c(254,255))))enc<-"UTF-16BE"
  else enc <- .io_or(encoding,"UTF-8")
  lines <- readLines(path,warn=FALSE,encoding=enc)
  lines <- sub("^\ufeff","",lines)
  if(skip_header>0L)lines<-lines[-seq_len(min(skip_header,length(lines)))]
  begin <- grep("begin spectral data",lines,ignore.case=TRUE)
  if(length(begin)) {
    end <- grep("end spectral data",lines,ignore.case=TRUE); end<-if(length(end))end[1L]-1L else length(lines)
    lines<-lines[seq.int(begin[1L]+1L,end)]
  }
  if(!is.null(comments)&&nzchar(comments)) {
    trim <- trimws(lines)
    keep <- !startsWith(trim,comments)|grepl("^#\\s*(wave|x|y|z|intensity|time)\\b",trim,ignore.case=TRUE,perl=TRUE)
    lines<-lines[keep]
  }
  lines<-lines[nzchar(trimws(lines))]
  delim <- .io_or(delimiter,default_delimiter)
  if(is.null(delim)&&length(lines)) {
    cand <- c("\t",";",",","|")
    scores <- vapply(cand,function(d) {
      counts<-lengths(strsplit(head(lines,400L),d,fixed=TRUE))-1L
      tb<-table(counts);if(!length(tb))return(0)
      tb<-tb[names(tb)!="0"];if(length(tb))max(tb)else 0
    },numeric(1))
    if(max(scores)>0)delim<-cand[which.max(scores)]
    if(identical(delim,",")&&all(grepl("^\\s*[+-]?[0-9]+,[0-9]+(\\s+[+-]?[0-9]+,[0-9]+)*\\s*$",tail(lines,20L))))delim<-NULL
  }
  rows <- lapply(lines,function(l) {
    if(is.null(delim))strsplit(trimws(l),"[[:space:]]+")[[1L]]
    else if(grepl('"',l,fixed=TRUE)) as.character(utils::read.table(text=l,sep=delim,header=FALSE,quote='"',comment.char="",fill=TRUE,colClasses="character",na.strings=NULL)[1,])
    else strsplit(paste0(l,delim),delim,fixed=TRUE)[[1L]]
  })
  comma <- identical(decimal,",")||(is.null(decimal)&&!identical(delim,",")&&any(grepl("^[+-]?[0-9]+,[0-9]+([eE][+-]?[0-9]+)?$",unlist(tail(rows,50L)))))
  if(comma)rows<-lapply(rows,function(r){i<-grepl("^[+-]?[0-9]+,[0-9]+([eE][+-]?[0-9]+)?$",r);r[i]<-sub(",",".",r[i],fixed=TRUE);r})
  grid<-.io_grid(rows)
  .io_table(grid$cols,grid$header,layout)
}
.io_long <- function(rec) {
  rows <- .io_rows(rec$arrays[[1L]]); lens<-lengths(rows); n<-length(rows)
  cols <- list(spectrum=rep(seq_len(n),lens))
  if(!is.null(rec$names))cols$name<-rep(rec$names,lens)
  if(!is.null(rec$positions))for(j in seq_len(ncol(rec$positions)))cols[[rec$position_names[j]]]<-rep(rec$positions[,j],lens)
  cols$wavenumber <- if(is.null(rec$wavenumbers))unlist(lapply(lens,function(l)seq_len(l)-1L))else if(is.list(rec$wavenumbers))unlist(rec$wavenumbers)else rep(rec$wavenumbers,n)
  for(k in names(rec$arrays))cols[[k]]<-unlist(.io_rows(rec$arrays[[k]]),use.names=FALSE)
  as.data.frame(cols,check.names=FALSE,stringsAsFactors=FALSE)
}
.io_wide <- function(rec) {
  if(length(rec$arrays)!=1L||is.list(rec$arrays[[1L]]))stop("wide layout requires one rectangular array",call.=FALSE)
  a<-rec$arrays[[1L]]
  out<-data.frame(wavenumber=.io_or(rec$wavenumbers,seq_len(ncol(a))-1L),t(a),check.names=FALSE)
  names(out)<-c(if(is.null(rec$wavenumbers))"index"else"wavenumber",make.unique(.io_or(rec$names,paste0("spectrum_",seq_len(nrow(a))))))
  out
}
.io_choose <- function(rec,layout=NULL,wide=TRUE) {
  if(is.null(layout))layout<-if(wide&&length(rec$arrays)==1L&&!is.list(rec$arrays[[1L]])&&is.null(rec$positions)&&nrow(rec$arrays[[1L]])<ncol(rec$arrays[[1L]]))"wide"else"long"
  if(layout=="wide").io_wide(rec)else if(layout=="long").io_long(rec)else stop("Writing layout must be wide or long",call.=FALSE)
}
.io_write_text <- function(path,rec,delimiter=",",decimal=".",layout=NULL,...) {
  cols<-.io_choose(rec,layout)
  old<-options(digits=17L);on.exit(options(old))
  utils::write.table(cols,path,sep=delimiter,dec=decimal,quote=TRUE,row.names=FALSE,col.names=TRUE,na="NaN",fileEncoding="UTF-8")
}
.io_read_excel <- function(path,sheet=NULL,layout="auto",...) {
  .io_need("readxl")
  sheets<-readxl::excel_sheets(path); sheets<-setdiff(sheets,"radar_metadata")
  if(!is.null(sheet)&&!identical(sheet,"all"))sheets<-if(is.numeric(sheet))sheets[sheet]else sheet
  # RADAR component workbooks use one sheet per array; select all automatically.
  if(is.null(sheet)&&!all(c("spectra",.io_components)%in%sheets))sheets<-sheets[1L]
  records<-lapply(sheets,function(s){a<-as.data.frame(readxl::read_excel(path,sheet=s,col_names=FALSE,.name_repair="minimal"));g<-.io_grid(lapply(seq_len(nrow(a)),function(i)as.character(a[i,])));.io_table(g$cols,g$header,layout)})
  rec<-records[[1L]]
  if(length(records)>1L)rec$arrays<-setNames(lapply(records,`[[`,"intensities"),sheets)
  if("radar_metadata"%in%readxl::excel_sheets(path)) {
    m<-readxl::read_excel(path,sheet="radar_metadata",col_names=FALSE,.name_repair="minimal")
    if(nrow(m))rec<-.io_apply_meta(rec,as.character(m[[1L]][nrow(m)]))
  }
  .io_record(rec)
}
.io_write_excel <- function(path,rec,layout=NULL,...) {
  .io_need("openxlsx")
  if(tolower(tools::file_ext(path))=="xls")stop("Writing .xls is unsupported; use .xlsx",call.=FALSE)
  if(identical(layout,"long")||is.list(rec$arrays[[1L]])||!is.null(rec$positions))sheets<-list(spectra=.io_long(rec))
  else if(length(rec$arrays)==1L)sheets<-list(spectra=.io_wide(rec))
  else sheets<-lapply(rec$arrays,function(a){r<-rec;r$arrays<-list(intensity=a);.io_wide(r)})
  sheets$radar_metadata<-data.frame(json=.io_meta_json(rec),check.names=FALSE)
  openxlsx::write.xlsx(sheets,path,overwrite=TRUE,colNames=TRUE)
}
.io_read_sqlite <- function(path,table=NULL,query=NULL,layout="auto",...) {
  .io_need("DBI");.io_need("RSQLite")
  con<-DBI::dbConnect(RSQLite::SQLite(),path);on.exit(DBI::dbDisconnect(con))
  tabs<-DBI::dbListTables(con)
  if(is.null(table))table<-if("spectra"%in%tabs)"spectra"else setdiff(tabs,c("radar_metadata","_radar_metadata"))[1L]
  cols<-if(!is.null(query))DBI::dbGetQuery(con,query)else DBI::dbReadTable(con,table)
  rec<-.io_table(cols,names(cols),layout)
  mt<-intersect(c("radar_metadata","_radar_metadata"),tabs)
  if(length(mt)){
    m<-DBI::dbReadTable(con,mt[1L])
    if(nrow(m)) {
      if(all(c("key","value")%in%names(m))) {
        .io_need("jsonlite")
        meta<-setNames(lapply(m$value,function(v)tryCatch(jsonlite::fromJSON(v,simplifyVector=FALSE),error=function(e)v)),m$key)
        rec<-.io_apply_meta(rec,meta)
      } else rec<-.io_apply_meta(rec,as.character(m[[ncol(m)]][1L]))
    }
  }
  rec
}
.io_write_sqlite <- function(path,rec,layout="long",table="spectra",if_exists="replace",...) {
  .io_need("DBI");.io_need("RSQLite")
  if(!if_exists%in%c("replace","append","fail"))stop("if_exists must be replace, append or fail",call.=FALSE)
  con<-DBI::dbConnect(RSQLite::SQLite(),path);on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con,table,.io_choose(rec,layout,FALSE),overwrite=if_exists=="replace",append=if_exists=="append")
  .io_need("jsonlite");m<-.io_meta(rec)
  metadata<-data.frame(key=names(m),value=vapply(m,function(v)as.character(jsonlite::toJSON(v,auto_unbox=TRUE,null="null",digits=NA)),character(1)))
  DBI::dbWriteTable(con,"radar_metadata",metadata,overwrite=TRUE)
}
.io_read_arrow <- function(path,feather=FALSE,layout="auto",...) {
  .io_need("arrow")
  tab<-if(feather)arrow::read_ipc_file(path,as_data_frame=FALSE)else arrow::read_parquet(path,as_data_frame=FALSE)
  rec<-.io_table(as.data.frame(tab),tab$ColumnNames,layout)
  js<-tab$schema$metadata$radar
  if(!is.null(js))rec<-.io_apply_meta(rec,js)
  rec
}
.io_write_arrow <- function(path,rec,feather=FALSE,layout=NULL,...) {
  .io_need("arrow")
  tab<-arrow::Table$create(.io_choose(rec,layout,FALSE))
  tab<-tab$ReplaceSchemaMetadata(c(tab$schema$metadata,list(radar=.io_meta_json(rec))))
  if(feather)arrow::write_ipc_file(tab,path)else arrow::write_parquet(tab,path)
}
