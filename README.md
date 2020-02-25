# rCOGS

This package contains routines to prioritise putative causal genes and tissue contexts using promoter-capture Hi-C data. COGS stands for __C__ apture Hi-C __O__ mnibus __G__ ene __S__ core and works best with pcHi-C datasets processed using [CHiCAGO](https://bioconductor.org/packages/release/bioc/html/Chicago.html). A full description of the approach can be found in [Javierre et al.](http://dx.doi.org/10.1016/j.cell.2016.09.037) and [Burren et al.](http://dx.doi.org/10.1186/s13059-017-1285-0).

# Installation

```
library(devtools)
install_github('ollyburren/rCOGS',build_vignettes=TRUE)
````

# Documentation

Please see online documentation and vignettes for common use cases [here](https://ollyburren.github.io/rCOGS).

A quickstart guide for performing your own analyses is available [here](https://ollyburren.github.io/rCOGS/articles/Quickstart.html)
