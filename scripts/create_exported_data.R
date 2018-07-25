library(devtools)
load_all('/Users/oliver/git/rCOGS')


## load set of summary statistics and compute ppi

ld.gr <- load_ld_regions('/Users/oliver/DATA/rCOGS/hapmap_recomb.bed')
maf.DT <- load_ref_maf('/Users/oliver/DATA/rCOGS/uk10k_0.001.tab')
sle <- load_gwas('/Users/oliver/DATA/rCOGS/SLE.tab',maf.DT,ld.gr,n.cases=4036,n.controls=6959)
## remove MHC

mhc.idx <- which(sle$chr==6 & between(sle$pos,25e6,35e6))
sle <- sle[-mhc.idx,]


regions <- make_pchic('/Users/oliver/DATA/rCOGS/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab',biotype.filter='protein_coding')
all.genes <- lapply(regions,function(g) unique(g$ensg)) %>% do.call('c',.) %>% unique
regions[['VProm']] <- make_vprom('/Users/oliver/DATA/rCOGS/Digest_Human_HindIII.tab','/Users/oliver/DATA/rCOGS/HindIII_baits_e75.tab',all.genes)
csnps <- make_csnps('/Users/oliver/DATA/rCOGS/cSNPs_w_ENSG.e75.tab')


## next compute indicator matrix but base this on fragments rather than genes
digest <- load_digest('/Users/oliver/DATA/rCOGS/Digest_Human_HindIII.tab')

fs.overall.scores <- compute_cogs(sle,csnps,digest,regions)

## compare these score with those from previous incarnations.

old <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/out/geneScore/SLE.pmi.tab')[biotype=='protein_coding',.(ensg,all_gene_score)]

M <- merge(old,fs.overall.score,by.x='ensg',by.y='ensg',all=TRUE)
setnames(M,c('ensg','old','new'))
library(cowplot)

ggplot(M,aes(old,new)) + geom_point()
