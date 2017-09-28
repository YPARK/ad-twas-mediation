
options(stringsAsFators = FALSE)
source('util.R')
source('figure.util.R')
library(dplyr)
library(ashr)

dir.create('figures', recursive = TRUE, showWarnings = FALSE)

read.mediation <- function(in.hdr, in.cols) {
    .files <- in.hdr %&&% 1:22 %&&% '.mediation.gz'
    .dat <- lapply(.files, read.table, col.names = in.cols)
    return(do.call(rbind, .dat))
}

expr.cols <- c('ensg', 'theta', 'theta.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.1', 'ld.2', 'n.snps')

data.direct.tab <- read.mediation('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_', expr.cols)
data.marginal.tab <- read.mediation('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_', expr.cols)

## bootstrap p-value distribution
plt.hist <- function(.df) {
    plt <- ggplot(.df, aes(x=p.val, group = null, fill = null)) + theme_bw() +
        geom_histogram(position = 'dodge', bins = 20, color = 'gray40') +
            geom_hline(yintercept = nrow(.df) / 2 / 20, lty = 2)

    plt <- plt +
        scale_fill_manual(values = c('#AAAAFF', 'gray80'))

    plt <- plt +
        theme(legend.title = element_blank(),
              legend.position = c(.9, .9),
              legend.justification = c(1, 1))
}

.df <- rbind(data.direct.tab %>% select(emp.p) %>% mutate(null = 'direct'),
             data.marginal.tab %>% select(emp.p) %>% mutate(null = 'marginal')) %>%
    mutate(p.val = 10^(-emp.p))

out.file <- 'figures/bootstrap_gene.pdf'
pdf(file = out.file, useDingbats = FALSE, width = 6, height = 4)
plt.hist(.df) + xlab('Empirical p-value')
dev.off()

.df <- rbind(data.direct.tab %>% select(norm.p) %>% mutate(null = 'direct'),
             data.marginal.tab %>% select(norm.p) %>% mutate(null = 'marginal')) %>%
    mutate(p.val = 10^(-norm.p))

out.file <- 'figures/bootstrap_gene_approx.pdf'
pdf(file = out.file, useDingbats = FALSE, width = 6, height = 4)
plt.hist(.df) + xlab('Approximate p-value')
dev.off()

## PIP cutoff vs bootstrap p-value
plt.pip.cutoff <- function(.df) {
    plt.aes <- aes(x = pip, y = mean.p,
                   ymin = pmax(mean.p - se.p, 0),
                   ymax = pmin(mean.p + se.p, 1))

    p1 <- ggplot(.df, plt.aes) + theme_bw() +
        geom_ribbon(fill = '#CCFFCC') + geom_point(size = 2) + geom_line(color = 'gray50')

    p1 <- p1 + xlab('Log-odds PIP') + ylab('Bootstrap P-value')
    return(p1)
}

pip.marginal.tab <- data.marginal.tab %>%
    mutate(pip = pmin(pmax(round(lodds * 5)/5, -5), 3)) %>%
    group_by(pip) %>%
    summarize(mean.p = mean(10^(-emp.p)), se.p = sd(10^(-emp.p)))

pip.direct.tab <- data.direct.tab %>%
    mutate(pip = pmin(pmax(round(lodds * 5)/5, -5), 3)) %>%
    group_by(pip) %>%
    summarize(mean.p = mean(10^(-emp.p)), se.p = sd(10^(-emp.p)))

out.file <- 'figures/bootstrap_gene_cutoff_marginal.pdf'
pdf(file = out.file, useDingbats = FALSE, width = 4, height = 4)
print(plt.pip.cutoff(pip.marginal.tab) + ggtitle('Null = marginal effect model'))
dev.off()

out.file <- 'figures/bootstrap_gene_cutoff_direct.pdf'
pdf(file = out.file, useDingbats = FALSE, width = 4, height = 4)
print(plt.pip.cutoff(pip.direct.tab) + ggtitle('Null = direct effect model'))
dev.off()

## PIP cutoff vs ASH q-value
ash.direct <- ash(betahat = data.direct.tab$theta, sebetahat = data.direct.tab$theta.se)
ash.marginal <- ash(betahat = data.marginal.tab$theta, sebetahat = data.marginal.tab$theta.se)

fdr.df <- cbind(data.direct.tab %>% select(lodds),
                       ash.direct$result %>% select(qvalue)) %>%
    mutate(pip = round(lodds * 5)/5) %>%
    group_by(pip) %>%
    summarize(fdr.max = max(qvalue), fdr.min = min(qvalue))

plt.pip.ash <- function(.df) {
    plt.aes <- aes(x = pip, ymax = fdr.max, ymin = fdr.min)
    p1 <- ggplot(.df, plt.aes) + theme_bw() +
        geom_ribbon(fill = '#CCFFCC')

    p1 <- p1 + geom_line(aes(y = fdr.max)) +
        geom_line(aes(y = fdr.min), lty = 2) +
        xlab('Log-odds PIP') + ylab('FDR')
            
    p1 <- p1 + geom_hline(yintercept = 0.1, color = 'red', lty = 2) +
        scale_x_continuous(breaks = seq(-5, 5)) +
            scale_y_continuous(breaks = seq(0, 1, .05))

    return(p1)
}

out.file <- 'figures/fdr_gene_cutoff.pdf'
pdf(file = out.file, useDingbats = FALSE, width = 4, height = 4)
print(plt.pip.ash(fdr.df))
dev.off()
