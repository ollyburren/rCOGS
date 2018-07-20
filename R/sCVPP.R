#' helper function to sum logs without loss of precision
#' \code{logsum} sums logs without loss of precision
#'
#' @param x a vector of logs to sum
#' @return a scalar

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

#' helper function to compute shrinkage factor for a quantitative trait
#' @param f a vector of length $i$ of minor allele frequencies in controls matching the order of p
#' @param N a scalar or vector (of length $i$) of the total number of samples included in association test
#' @return a vector of shrinkage factors


Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

#' helper function to compute shrinkage factor for a case/control trait
#' @param f a vector of length $i$ of minor allele frequencies in controls matching the order of p
#' @param N a scalar or vector (of length $i$) of the total number of samples included in association test
#' @param s a scalar or vector (of length $i$) of the proportion of cases for each test.
#' @return a vector of shrinkage factors

## compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}


#' This function computes posterior probabilities for a variant to be causal using Wakefields asymptotic Bayes Factors.
#' @param p vector of p.values for $i$ univariate tests of association for the trait of interest with a set of variants
#' @param f a vector of length $i$ of minor allele frequencies in controls matching the order of p
#' @param type a scalar taking the  of the proportion of cases for each test.
#' @param N a scalar or vector (of length $i$) of the total number of samples included in association test
#' @param s a scalar or vector (of length $i$) of the proportion of cases for each test.
#' @param pi_i a scalar of the prior probability that a given variant is causal (default=1e-4)
#' @return a vector of prior probabilities.
#' @export


## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f,type=c('QUANT','CC'), N, s,pi_i=1e-4) {
  if(type=="QUANT") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else if (type=='CC') {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }else{
    stop(sprintf("COGSR:approx.bf.p invalid study type expecting 'QUANT' or 'CC' got '%s'",type))
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- sd.prior^2 / (sd.prior^2 + V)
  lABF = 0.5 * (log(1-r) + (r * z^2))
  sBF <- logsum(lABF + log(pi_i))
  ppi<- exp(lABF + log(pi_i))/(exp(sBF) + 1)
  return(ppi)
}
