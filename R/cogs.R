#' This function computes COGS scores for a given trait across a set of pchic regions
#' @param gwas - data.table obtained from load_gwas
#' @param csnps = GRanges object obtained from make_csnps
#' @param digest GRanges object obtained from load_digest
#' @param regions list object obtained from make_pchic
#' @param feature.names - a character vector of regions to include these must be present in pcHi-C file
#' @return data.table of results
#' @export

compute_cogs <- function(gwas,csnps,digest,regions,feature.names){
  gwas.gr <- with(gwas,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L)))
  ## annotate coding SNPs - we remove these and add back in for gene of interest
  ol.csnps <- findOverlaps(csnps,gwas.gr) %>% as.matrix
  gwas[ol.csnps[,2],csnp:=csnps[ol.csnps[,1],]$ensg]
  csnps.DT <- gwas[,list(csnp.ppi=sum(ppi)),by=c('csnp','ld')][!is.na(csnp),][,.(ensg=csnp,ld=ld,frag.ppi=csnp.ppi)]
  ol <- findOverlaps(digest,gwas.gr) %>% as.matrix
  ol.DT <- data.table(gwas[ol[,2],.(csnp,ppi,ld)],fragid=digest[ol[,1],]$fragid)
  #next we remove cSNPs and compute sum for each fragment
  ol.DT <- ol.DT[is.na(csnp),][,list(frag.ppi=sum(ppi)),by=c('ld','fragid')]
  ## to compute scores we want for each lu[['ensg']][['tissue']] is a vector of fragids
  all <- lapply(names(regions),function(ct){
    tDT <- regions[[ct]]
    tDT[,.(ensg,fragid,tissue=ct)]
  }) %>% rbindlist
  mDT <- merge(all,ol.DT,by.x='fragid',by.y='fragid')
  mDT[,tissue:=factor(tissue)]
  ## we can compute any score by selecting one or more tissues
  if(missing(feature.names)){
      feature.names <- c(names(regions),'coding_snp')
      message(sprintf("rCOGS:compute_cogs feature.names argument missing computing overall score using %s",paste(feature.names,collapse=",\n",sep=",\n")))
  }
  #feature.names <- c(names(regions),'coding_snp')
  fs.DT <- mDT[tissue %in% feature.names,]
  ## remove duplicate fragments by gene so don't count twice
  fs.DT <- fs.DT[!duplicated(fs.DT[,.(fragid,ensg)]),][,.(ensg,ld,frag.ppi,fragid)]
  ## check to see if coding snp is selected if so add
  #if(sum(feature.names=='coding_snp')>0)
  if(any(feature.names=='coding_snp'))
    fs.DT <- rbind(fs.DT,csnps.DT[,fragid:=0])
  fs.ld.score <- fs.DT[,list(ld.score=sum(frag.ppi)),by=c('ld','ensg')]
  fs.ld.score[,list(cogs=1-prod(1-ld.score)),by='ensg']
}
