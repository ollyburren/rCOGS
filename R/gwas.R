#' @import data.table


#library(GenomicRanges)
#library(data.table)
#library(magrittr)




MAF_COLS <- c('chr','pos','maf')
LD_COLS <- c('chr','start','end')


#' This function fixes chromosomes so they are numeric
#' @param chr a vector of strings representing chromosomes
#' @param f a scalar representing a file name (so user can understand which file is problematic)
#' @return a vector of reformatted chromosomes

fix_chr <- function(chr,f){
  ct <- gsub("^chr","",chr)
  ct[ct=='X']=23
  ct[ct=='Y']=24
  ct[ct=='MT']=25
  if(grepl("^[0-9][0-9]*$",ct) %>% sum != length(ct))
    stop(sprintf("File %s contains non standard chromosomes",f))
  return(as.numeric(ct))
}

#' This function checks and loads a reference set of maf expects three columns tab delim file with colnames  chr,pos and maf
#' @param f character scalar representing filename containing maf data expects three columns wtih names chr, pos and maf.
#' @param min.maf scalar numeric representing minor allele threshold cut off.
#' @return data.table object of minor allele freq.
#' @export


load_ref_maf <- function(f,min.maf=0.05){
  DT <- fread(f)
  if(!identical(names(DT),MAF_COLS)){
    stop(sprintf("COGSR:load_ref_maf file format error, expected columns ''%s' got '%s' ",
      paste(MAF_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  ## check for allele freq rather than maf and fix accordingly
  DT[,maf:=ifelse(maf>0.5,1-maf,maf)]
  message(sprintf("Adding unique identifier and filtering on MAF > %f",min.maf))
  DT <- DT[maf>min.maf][,pid:=fix_chr(chr,f) %>% paste(.,pos,sep=':')]
  setkey(DT,pid)
  return(DT)
}

#' This function checks and loads and checks a set of non overlapping  LD regions expects three column tab delim file with colnames chr,start,end
#' @param f character scalar representing filename containing maf data expects three columns wtih names chr, start and end.
#' @return GenomicRanges object of ld blocks/regions
#' @export

load_ld_regions <- function(f){
  DT <- fread(f)
  if(!identical(names(DT),LD_COLS)){
    stop(sprintf("COGSR:load_ld_regions file format error, expected columns ''%s' got '%s'",
      paste(LD_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  #DT[,chr:=fix_chr(chr,f)]
  DT <- DT[order(chr,start),][,ldid:=1:.N]
  gr <- with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ldid=ldid))
  ## do some checks to see that these are non overlapping regions
  if(length(reduce(gr,min.gapwidth=0))!=length(gr)){
    stop(sprintf("OGSR:load_ld_regions file format error, ld regions appear to overlap please check"))
  }
  gr
}

#' This function loads and checks set of univariate P values for a given trait and adds posterior probabilities for each SNP to be causal
#' using Wakefield's aBF
#' @param f character scalar representing GWAS file to load, tab delimited with three columns chr,pos and p
#' @param maf.DT a data.table of reference MAFs created by load_ref_maf
#' @param ld.gr a GenomicRanges object of LD regions created by load_ld_regions
#' @param n.cases a scalar for the number of GWAS cases - this should equal sample size in the case of a quantitative trait
#' @param n.controls a scalar for the number of GWAS controls - this should be omitted in the case of a quantitative trait
#' @param plink a boolean as to whether input format is PLINK rather than the cannonical type described above.
#' @return a data.table of annotated GWAS variants suitable for COGS analysis
#' @export

load_gwas <- function(f,maf.DT,ld.gr,n.cases,n.controls,plink=FALSE){
  DT <- fread(f)
  if(plink){
    DT <- DT[,.(CHR,BP,P)]
    setnames(DT,c('chr','pos','p'))
  }
  message("Adding unique identifier")
  DT[,pid:=fix_chr(chr,f) %>% paste(.,pos,sep=':')]
  setkey(DT,pid)
  ## add maf and filter
  M.DT <- merge(DT,maf.DT[,.(pid,maf)],by.x='pid',by.y='pid')
  message(sprintf("COGSR:load_gwas added MAF before %d variants after %d",nrow(DT),nrow(M.DT)))
  ## add ld block information
  M.gr <- with(M.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),pid=pid))
  ol <- findOverlaps(M.gr,ld.gr) %>% as.matrix()
  M.DT[ol[,1],ld:=ol[,2]]
  ld.miss <- sum(is.na(M.DT$ld))
  if(ld.miss>0){
    message(sprintf("COGSR:load_gwas %d input variants did not overlap an LD block these won't be included in analysis",ld.miss))
    M.DT <- M.DT[!is.na(ld),]
  }
  if(any(duplicated(DT$M.DT)))
      stop(sprintf("COGSR:load_gwas found duplicate variants by position in %s pleas fix",f))
  if(missing(n.controls)){
    message("Control set sample size not set assuming a quantitative trait")
    M.DT[,ppi:=approx.bf.p(p,maf,'QUANT',n.cases),by=ld]
  }else{
    study.size <- n.cases + n.controls
    prop <- n.cases/study.size
    M.DT[,ppi:=approx.bf.p(p,maf,'CC',study.size,prop),by=ld]
  }
  ## check to see if any ppi are NA as this indicates an error
  na.ppi <- sum(is.na(M.DT$ppi)) %>% sum
  if(na.ppi>0)
      message(sprintf("COGSR:load_gwas %d input variants have an invalid ppi please check that ld block",na.ppi))
  M.DT
}
