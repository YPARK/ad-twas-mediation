
options(stringsAsFators = FALSE)
source('util.R')
source('figure.util.R')
library(dplyr)

read.mediation <- function(in.hdr, in.cols) {
    .files <- in.hdr %&&% 1:22 %&&% '.mediation.gz'
    .dat <- lapply(.files, read.table, col.names = in.cols)
    return(do.call(rbind, .dat))
}

expr.cols <- c('ensg', 'theta', 'thet.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.1', 'ld.2', 'n.snps')

data.direct.tab <- read.mediation('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_', expr.cols)
data.marginal.tab <- read.mediation('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_', expr.cols)


## take.hist <- function

ggplot(data.direct.tab %>% filter(lodds < 0)) +
    geom_histogram(aes(x=10^(-emp.p)))

ggplot(data.marginal.tab %>% filter(lodds < 0)) +
    geom_histogram(aes(x=10^(-emp.p)))

