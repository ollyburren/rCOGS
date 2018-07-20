library(devtools)
load_all('/Users/oliver/git/rCOGS')


## load set of summary statistics and compute ppi

ld.gr <- load_ld_regions('/Users/oliver/DATA/rCOGS/hapmap_recomb.bed')
maf.DT <- load_ref_maf('/Users/oliver/DATA/rCOGS/uk10k_0.001.tab')
sle <- load_gwas('/Users/oliver/DATA/rCOGS/SLE.tab',maf.DT,ld.gr,n.cases=4036,n.controls=6959)
regions <- make_pchic('/Users/oliver/DATA/rCOGS/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab',biotype.filter='protein_coding')
all.genes <- lapply(regions,function(g) unique(g$ensg)) %>% do.call('c',.) %>% unique
regions[['VProm']] <- make_vprom('/Users/oliver/DATA/rCOGS/Digest_Human_HindIII.tab','/Users/oliver/DATA/rCOGS/HindIII_baits_e75.tab',all.genes)
csnps <- make_csnps('/Users/oliver/DATA/rCOGS/cSNPs_w_ENSG.e75.tab')


## next compute indicator matrix but base this on fragments rather than genes
digest <- load_digest('/Users/oliver/DATA/rCOGS/Digest_Human_HindIII.tab')
gwas.gr <- with(sle,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L)))
## annotate coding SNPs - we remove these and add back in for gene of interest
ol.csnps <- findOverlaps(csnps,gwas.gr) %>% as.matrix
sle[ol.csnps[,2],csnp:=csnps[ol.csnps[,1],]$ensg]
csnps.DT <- sle[,list(csnp.ppi=sum(ppi)),by=c('csnp','ld')][!is.na(csnp),][,.(ensg=csnp,ld=ld,frag.ppi=csnp.ppi)]


ol <- findOverlaps(digest,gwas.gr) %>% as.matrix
ol.DT <- data.table(sle[ol[,2],.(csnp,ppi,ld)],fragid=digest[ol[,1],]$fragid)
#next we remove cSNPs and compute sum for each fragment
ol.DT <- ol.DT[is.na(csnp),][,list(frag.ppi=sum(ppi)),by=c('ld','fragid')]

## to compute scores we want for each lu[['ensg']][['tissue']] is a vector of fragids

all <- lapply(names(regions),function(ct){
  tDT <- regions[[ct]]
  tDT[,.(ensg,fragid,tissue=ct)]
}) %>% rbindlist

mDT <- merge(all,ol.DT,by.x='fragid',by.y='fragid')
#mDT <- rbind(mDT,csnps)
mDT[,tissue:=factor(tissue)]

## we can compute any score by selecting one or more tissues
feature.set <- c(names(regions),'coding_snp')
fs.DT <- mDT[tissue %in% feature.set,]
## remove duplicate fragments by gene so don't count twice
fs.DT <- fs.DT[!duplicated(fs.DT[,.(fragid,ensg)]),][,.(ensg,ld,frag.ppi)]
## check to see if coding snp is selected if so add
if(sum(feature.set=='coding_snp')>0)
  fs.DT <- rbind(fs.DT,csnps.DT)
fs.ld.score <- fs.DT[,list(ld.score=sum(frag.ppi)),by=c('ld','ensg')]
fs.overall.score <- fs.ld.score[,list(cogs=1-prod(1-ld.score)),by='ensg']
