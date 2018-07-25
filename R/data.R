#' ld.gr
#' HapMap approx 1cM LD blocks for chromosome 22.
#'
#' @format GRanges object with 69 ranges and 1 metadata column
"ld.gr"

#' maf.DT
#' Reference minor allele frequencies from UK10K for chromosome 22.
#'
#' @format data.table with 74356 rows and 4 variables
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{pos}{position GRCh37}
#'   \item{maf}{minor allele frequency}
#'   \item{pid}{unique id obtained from chr and position}
#' }
"maf.DT"

#' t1d.DT
#' GWAS data from Cooper et al. for Type 1 diabetes for chromosome 22.
#'
#' @format data.table with 68300 rows and 7 variables
#' \describe{
#'   \item{pid}{unique id obtained from chr and position}
#'   \item{chr}{chromosome}
#'   \item{pos}{position GRCh37}
#'   \item{p}{univariate p value}
#'   \item{maf}{minor allele freq from maf.DT}
#'   \item{ld}{ld block from ld.gr}
#'   \item{ppi}{posterior probability for variant to be causal}
#' }
"t1d.DT"

#' feature.sets
#' Feature sets from Javierre et al pcHi-C dataset for chromosome 22.
#'
#' @format list of 18 data.table objects
"feature.sets"


#' csnps.gr
#' Coding SNPs for chromosome 22.
#'
#' @format GRanges object with 53002 ranges and 2 metadata columns.
"csnps.gr"


#' digest.gr
#' HindIII restriction digest for chromosome 22.
#'
#' @format GRanges object with 7582 ranges and 1 metadata column.
"digest.gr"
