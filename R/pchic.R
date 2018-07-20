## routine to load promoter pcHi-C DATA

library(data.table)
library(magrittr)
library(GenomicRanges)

PCHIC_COLS <- c("ensg","name","biotype","strand","baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist")
DIGEST_COLS <- c('chr','start','end','fragid')
DESIGN_COLS <- c('fragid','ensg')
CSNP_COLS <- c('chr','pos','ensg')

#' This loads a restriction digest file with columns chr,start,end and fragid
#' @param f a scalar representing a file name of digest file
#' @return GRanges object representing restriction digest
#' @export

load_digest <- function(f){
  DT <- fread(f)
  if(!identical(names(DT),DIGEST_COLS)){
    stop(sprintf("COGSR:load_digest file format error, expected columns ''%s' got '%s'",
      paste(DIGEST_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT[,chr:=fix_chr(chr,f)]
  DT <- DT[order(chr,start),][chr!=25,]
  gr <- with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),fragid=fragid))
  if(length(reduce(gr,min.gapwidth=0))!=length(gr)){
    stop(sprintf("OGSR:load_ld_regions file format error, ld regions appear to overlap please check"))
  }
  gr
}

#' This loads a capture design file with columns fragid and ensg. fragid should match digest file.
#' @param f a scalar representing a file name of design file
#' @return data.table object representing capture design

load_capture_design <- function(f){
  DT <- fread(f)
  if(!identical(names(DT),DESIGN_COLS)){
    stop(sprintf("COGSR:load_capture_design file format error, expected columns ''%s' got '%s'",
      paste(DESIGN_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT
}

#' This function computes 'Virtual' promoter regions. These are baited fragments and the fragments either side.
#' @param f.digest a scalar representing a file name of digest file
#' @param f.design a scalar representing a file name of design file
#' @param keep.ensg a vector of ensembl id's to include. If ommited then all ensembl ids regardless of further filtering are included
#' @return data.table with two columns mapping ensg to virtual promoter fragment ids
#' @export

make_vprom <- function(f.digest,f.design,keep.ensg){
  frags.gr <- load_digest(f.digest)
  if(!missing(keep.ensg)){
    des.DT <- load_capture_design(f.design)[ensg %in% keep.ensg]
  }else{
    des.DT <- load_capture_design(f.design)
  }
  cap.gr <- frags.gr[frags.gr$fragid %in% des.DT$fragid,]
  lu.DT <- data.table(fragid=frags.gr$fragid,chr=seqnames(frags.gr) %>% as.character,start=start(frags.gr),end=end(frags.gr))
  ol<-findOverlaps(cap.gr,frags.gr,maxgap=1L) %>% as.matrix
  vprom.DT <- data.table(baitID=cap.gr[ol[,1],]$fragid,vpromID=frags.gr[ol[,2],]$fragid)
  vprom.DT<-merge(des.DT,vprom.DT,by.x='fragid',by.y='baitID')
  return(vprom.DT[,.(ensg,fragid=vpromID)])
  ## there is some sharing add duplicate rows
  # vprom.DT <- lu.DT[ol[,2],][,vprom:=ol[,1]]
  # vprom.DT <- vprom.DT[,list(start=min(start),end=max(end),count=.N),by=c('chr','vprom')][,fragid:=cap.gr[vprom,]$fragid]
  # if(any(vprom.DT$count!=3))
  #   stop(sprintf("COGSR:make_vprom error creating vproms cannot find adjacent fragments for some capture probes please check digest file %s and design file %s",f.digest,f.design))
  # vprom.DT <- merge(des.DT,vprom.DT,by.x='fragid',by.y='fragid')
  # with(vprom.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),fragid=fragid,ensg=ensg))
}

PCHIC_COLS <- c("ensg","name","biotype","strand","baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist")

#' This function computes loads in a set of promoter capture Hi-C files (pcHi-C). This is a tab delimited file that contains the following mandatory columns
#' \enumerate{
#' \item ensg - Ensembl gene id for gene
#' \item name - Gene name
#' \item biotype - Gene biotype from ensembl e.g. protein_coding
#' \item strand - Gene strand
#' \item baitChr - Bait (A captured restriction fragment) chromosome
#' \item baitStart - Bait fragment start
#' \item baitEnd - Bait fragment end
#' \item baitID - Bait fragment ID should match ID in digest file
#' \item baitName - Free text of what is captured - deprecated.
#' \item oeChr - 'Other end' (Promoter interacting restriction fragment) chromosome
#' \item oeStart - 'Other end' fragment start
#' \item oeEnd - 'Other end' fragment end
#' \item oeID - 'Other end' fragment ID should match ID in digest file
#' \item oeName - Free text of what is interacting - deprecated.
#' }
#' Additional columns contain CHiCAGO scores for one or more columns.
#' @param f a scalar representing file name of pcHi-C file
#' @param chic_thresh scalar representing CHiCAGO threshold for calling a significant interaction. Default is set to 5.
#' @param biotype.filter scalar/vector of characters representing a biotype filter. Default is set to 'protein_coding'
#' @param remove.trans.contacts a boolean on whether to remove interaction between chromosomes. Default is set to TRUE.
#' @return list of data.tables where each data.table has three columns ensg,fragid and chic.score
#' @export

make_pchic <- function(f,chic_thresh=5,biotype.filter='protein_coding',remove.trans.contacts=TRUE){
  DT <- fread(f)
  cnames <- names(DT)
  idx.ct.cols <- which(cnames %in% PCHIC_COLS)
  if(length(idx.ct.cols) != length(PCHIC_COLS))
    stop(sprintf("COGSR:load_pchic Missing required columns in pcHi-C file %s",f))
  if(!missing(biotype.filter)){
      message(sprintf("Filtering by biotypes %s",paste(biotype.filter,sep=",",collapse=',')))
      DT <- DT[biotype %in% biotype.filter,]
  }
  trans.idx <- which(DT$baitChr != DT$oeChr)
  if(remove.trans.contacts & length(trans.idx)>0){
    message(sprintf("Removing %d trans contacts",length(trans.idx)))
    DT <- DT[-trans.idx]
  }
  cell.types <- cnames[!cnames %in% PCHIC_COLS]
  DT <- DT[,c('baitChr','oeChr'):=list(fix_chr(baitChr),fix_chr(oeChr))]
  message(sprintf("Processing %d cell types",length(cell.types)))
  cts <- lapply(cell.types,function(ct){
    message(sprintf("Processing %s",ct))
    tDT <- DT[get(`ct`)>chic_thresh,.(ensg,fragid=oeID,chic.score=get(`ct`))]
    #tDT <- DT[get(`ct`)>chic_thresh,.(ensg,name,oeChr,oeStart,oeEnd,oeID,chic.score=get(`ct`))]
    #with(tDT,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),fragid=oeID,ensg=ensg))
  })
  names(cts) <- cell.types
  return(cts)
}


#' This function loads in a coding SNPs matched to a given ensembl id, columns are chr,pos and ensg
#' @param f a scalar representing file name of coding SNP file
#' @param keep.ensg a vector of ensembl id's to include. If ommited then all ensembl ids regardless of further filtering are included
#' @return GRanges object of coding SNPs
#' @export

make_csnps <- function(f,keep.ensg){
  if(!missing(keep.ensg)){
    DT <- fread(f)[ensg %in% keep.ensg]
  }else{
    DT <- fread(f)
  }
  if(!identical(names(DT),CSNP_COLS)){
    stop(sprintf("COGSR:load_capture_design file format error, expected columns ''%s' got '%s'",
      paste(CSNP_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT <- DT[order(chr,pos,ensg),]
  with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),fragid=0,ensg=ensg))
}
