
options(stringsAsFators = FALSE)
source('util.R')
source('figure.util.R')
library(dplyr)

out.file <- 'figures/bootstrap_gene.pdf'
dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

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

.df <- rbind(data.direct.tab %>% select(emp.p) %>% mutate(null = 'direct'),
             data.marginal.tab %>% select(emp.p) %>% mutate(null = 'marginal'))

plt <- ggplot(.df, aes(x=10^(-emp.p), group = null, fill = null)) + theme_bw() +
    geom_histogram(position = 'dodge', bins = 20, color = 'gray40') +
    geom_hline(yintercept = nrow(.df) / 2 / 20, lty = 2)

plt <- plt +
    xlab('Empirical p-value') +
    scale_fill_manual(values = c('#AAAAFF', 'gray80'))

plt <- plt +
    theme(legend.title = element_blank(),
          legend.position = c(.9, .9),
          legend.justification = c(1, 1))

pdf(file = out.file, useDingbats = FALSE, width = 6, height = 4)
print(plt)
dev.off()
