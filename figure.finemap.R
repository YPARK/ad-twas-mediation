argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

txn.stat.file <- argv[1]  # e.g., txn.stat.file <- 'finemap/IGAP_rosmap_eqtl_hs-lm_19.mediation.gz'
cpg.stat.file <- argv[2]  # e.g., cpg.stat.file <- 'finemap/IGAP_rosmap_mqtl_hs-lm_19.mediation.gz'
m2t.stat.file <- argv[3]  # e.g., m2t.stat.file <- 'm2t/IGAP_rosmap_hs-lm_19.mediation.gz'

out.hdr <- argv[4]

dir.create(dirname(out.hdr), recursive = TRUE)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ashr)
library(readr)
options(stringsAsFactors = FALSE)
source('figure.util.R')
source('util.R')

txn.col.names <- c('ensg', 'theta', 'theta.var', 'lodds',
                   'theta.2', 'theta.se.2', 'lodds.2',
                   'hgnc', 'tss', 'tes', 'strand',
                   'gwas.theta', 'gwas.z', 'chr', 'ld.start', 'ld.end',
                   'n.snps')

txn.stat <- read.table(txn.stat.file, col.names = txn.col.names) %>%
    arrange(tss) %>%
        mutate(loc = (tss+tes)/2/1e6, theta.se = sqrt(theta.var),
               gwas.lnp = pmin(-log10(2 * pnorm(abs(gwas.z), lower.tail = FALSE)), 15))

cpg.col.names <- c('cg', 'theta', 'theta.var', 'lodds',
                   'theta.2', 'theta.se.2', 'lodds.2',
                   'cg.loc', 'gwas.theta', 'gwas.z', 'chr', 'ld.start', 'ld.end',
                   'n.snps')

cpg.stat <- read.table(cpg.stat.file, col.names = cpg.col.names) %>%
    arrange(cg.loc) %>%
        mutate(loc = (cg.loc)/1e6, theta.se = sqrt(theta.var),
               gwas.lnp = pmin(-log10(2 * pnorm(abs(gwas.z), lower.tail = FALSE)), 15))

m2t.col.names <- c('ensg', 'cg', 'theta', 'theta.var', 'lodds',
                   'theta.2', 'theta.se.2', 'lodds.2',
                   'hgnc', 'tss', 'tes', 'strand',
                   'cg.loc', 'chr', 'ld.start', 'ld.end',
                   'n.snps')

m2t.stat <- read_tsv(m2t.stat.file, col_names = m2t.col.names)

txn.ash.out <- ash(txn.stat$theta, txn.stat$theta.se)
cpg.ash.out <- ash(cpg.stat$theta, cpg.stat$theta.se)

txn.stat <- cbind(txn.stat, qval = pmin(-log10(txn.ash.out$result$qvalue), 15))
cpg.stat <- cbind(cpg.stat, qval = pmin(-log10(cpg.ash.out$result$qvalue), 15))

fdr10.ensg <- txn.stat %>% filter(qval > 1) %>% select(ensg) %>% unique()
fdr10.cg <- cpg.stat %>% filter(qval > 1) %>% select(cg) %>% unique()

x.range <- range(cpg.stat$cg.loc, txn.stat$gene.loc)

plot.manhattan <- function(.df, font.size = 3) {

    .df$variable <- factor(.df$variable, c('mediation', 'gwas'), c('Mediation', 'GWAS'))
    .df.sig <- .df %>% filter(fdr > 1, variable == 'Mediation')
    
    ret <- ggplot() + theme_classic() + xlab('Genomic Location (Mb)') + ylab('Significance (-log10)') +
        geom_hline(yintercept = 2, color = '#FFCCCC') +
    geom_hline(yintercept = 8, color = 'gray60')

    ret <- ret +
        geom_point(data = .df %>% filter(variable == 'GWAS'), aes(x = loc, y = fdr),
                   size = 1, pch = 19, color = 'green')

    ret <- ret +
        geom_point(data = .df %>% filter(variable == 'Mediation'),
                   aes(x = loc, y = fdr, fill = sign(theta), size = 1/(1+exp(-lodds))),
                   pch = 21)
    ret <- ret +
        scale_fill_continuous(low = 'blue', high = 'yellow', guide = FALSE) +
            scale_size(range=c(0.1, 2), limits = c(0, 1), guide = FALSE)

    ret <- ret +
        geom_text_repel(data = .df.sig, aes(x = loc, y = fdr, label = label),
                        size = font.size)

    ret <- ret + ylim(c(0, 16))
                        
    return(ret)
}


.df <- cpg.stat %>%
    select(loc, theta, lodds, qval, gwas.lnp, cg) %>%
    rename(mediation = qval, gwas = gwas.lnp, label = cg) %>%
    melt(factorsAsStrings=TRUE, id.vars = c('loc','label', 'theta', 'lodds'),
         value.name = 'fdr') 

p0 <- plot.manhattan(.df, font.size = 2)


m2t.stat.linked <- m2t.stat %>% filter(lodds > 0)

.df <- cpg.stat %>% filter(cg %in% m2t.stat.linked$cg) %>%
    select(loc, theta, lodds, qval, gwas.lnp, cg) %>%
    rename(mediation = qval, gwas = gwas.lnp, label = cg) %>%
    melt(factorsAsStrings=TRUE, id.vars = c('loc','label', 'theta', 'lodds'),
         value.name = 'fdr') 

p1 <- plot.manhattan(.df)

m2t.stat.sig <- m2t.stat.linked %>% filter(cg %in% fdr10.cg[, 1]) %>%
    mutate(gene.loc = ifelse(strand == '+', tss, tes)) %>%
    arrange(cg.loc)

m2t.stat.sig <- m2t.stat.sig %>%
    left_join(m2t.stat.sig %>% select(cg) %>% unique() %>% mutate(cg.pos = 1:n()),
              by = 'cg') %>%
    left_join(m2t.stat.sig %>% select(ensg) %>% unique() %>% mutate(ensg.pos = 1:n()),
              by = 'ensg')

m2t.stat.sig <- m2t.stat.sig %>% mutate(cg.pos = cg.pos / max(cg.pos) * x.range[2]) %>%
    mutate(ensg.pos = ensg.pos / max(ensg.pos) * x.range[2])


p2.1 <-
    ggplot(m2t.stat.sig) +
    theme_void() +
    geom_segment(aes(x = cg.loc/1e6, xend = cg.pos/1e6, y = 1, yend = 0),
                 color = 'gray')

p2.2 <-
    ggplot(m2t.stat.sig) +
    theme_void() +
    geom_segment(aes(x = cg.pos/1e6, xend = ensg.pos/1e6, y = 1, yend = 0,
                     color = theta/sqrt(theta.var)),
                 arrow = arrow(length = unit(.25, 'lines'))) +
    scale_color_continuous(low = 'blue', high = 'red', guide = FALSE) +
    geom_point(aes(x = cg.pos/1e6, y = 1), size = 1) +
    geom_point(aes(x = ensg.pos/1e6, y = 0), size = 1)

p2.2 <- p2.2 +
    geom_text_repel(data = m2t.stat.sig %>% select(cg.pos, cg) %>% unique(),
                    aes(x = cg.pos/1e6, y = 1, label = cg), size = 3) +
    geom_text_repel(data = m2t.stat.sig %>% select(ensg.pos, hgnc) %>% unique(),
                    aes(x = ensg.pos/1e6, y = 0, label = hgnc), size = 3)

p2.3 <-
    ggplot(m2t.stat.sig) +
    theme_void() +
    geom_segment(aes(x = ensg.pos/1e6, xend = gene.loc/1e6, y = 1, yend = 0),
                 color = 'gray')


.df <- txn.stat %>% select(loc, qval, gwas.lnp, hgnc, theta, lodds) %>%
    rename(mediation = qval, gwas = gwas.lnp, label = hgnc) %>%
    melt(factorsAsStrings=TRUE, id.vars = c('loc','label','theta','lodds'),
         value.name = 'fdr')

p3 <- plot.manhattan(.df)

p1 <-
    p1 + geom_vline(data = m2t.stat.sig, aes(xintercept = cg.loc/1e6), color = 'gray', lty = 3)

p3 <-
    p3 + geom_vline(data = m2t.stat.sig, aes(xintercept = gene.loc/1e6), color = 'gray', lty = 3)

p.list <- lapply(list(p1, p2.1, p2.2, p2.3, p3),
                 function(xx) xx + scale_x_continuous(limits = x.range/1e6, expand = c(.01, .01)))


ww <- 1 + (x.range[2] - x.range[1]) / 1e6 / 10
hh <- 8
out.file <- out.hdr %&&% '-global.pdf'

pdf(file = out.file, useDingbats = FALSE, width = ww, height = hh)
grid.vcat(p.list, heights = c(3, .5, 1, .5, 3))
dev.off()

log.msg('%s\n\n', out.file)


out.file <- out.hdr %&&% '-global-cpg.pdf'

pdf(file = out.file, useDingbats = FALSE, width = ww, height = 3)
print(p0)
dev.off()

log.msg('%s\n\n', out.file)

################################################################
## show the best examples -- just connected ones
read.bedtools <- function(ld.bed.tab, tot.stat.file) {
    ## tot.stat.file <- 'stat/IGAP/data/hs-lm/19.eqtl_bed.gz'

    temp.dir <- system('mktemp -d', intern = TRUE, ignore.stderr = TRUE)
    bedtools.input <- temp.dir %&&% '/input.ucsc_bed.gz'
    bedtools.result <- temp.dir %&&% '/stat.ucsc_bed.gz'
    write.tab(ld.bed.tab, file = gzfile(bedtools.input))

    if(file.exists('/broad/software/scripts/useuse')) {
        bedtools.hdr <- 'if [ -f /broad/software/scripts/useuse ]; then source /broad/software/scripts/useuse > /dev/null; reuse -q BEDTools; fi; bedtools'
    } else {
        bedtools.hdr <- '/usr/local/bin/bedtools'
    }

    bedtools.cmd <- paste('intersect', '-a', tot.stat.file, '-b',
                          bedtools.input, ' -wa | gzip > ', bedtools.result)
    sum.cmd <- sprintf('%s %s', bedtools.hdr, bedtools.cmd)
    flag = system(sum.cmd)

    ret <- read_tsv(bedtools.result, col_names = FALSE)
}

ld.bed.tab <- m2t.stat.sig %>% select(chr, ld.start, ld.end) %>% unique()

stat.file <- 'stat/IGAP/data/hs-lm/' %&&% ld.bed.tab$chr[1] %&&% '.eqtl_bed.gz'
eqtl.stat <- read.bedtools(ld.bed.tab, stat.file)

colnames(eqtl.stat) <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
                         'a1', 'a2', 'qtl.theta', 'qtl.z',
                         'gwas.theta', 'gwas.se', 'gwas.z',
                         'ld', 'ensg', 'hgnc', 'tss', 'tes', 'strand')
                         
stat.file <- 'stat/IGAP/data/hs-lm/' %&&% ld.bed.tab$chr[1] %&&% '.mqtl_bed.gz'
mqtl.stat <- read.bedtools(ld.bed.tab, stat.file)

colnames(mqtl.stat) <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
                         'a1', 'a2', 'qtl.theta', 'qtl.z',
                         'gwas.theta', 'gwas.se', 'gwas.z',
                         'ld', 'cg', 'cg.loc')

## for each LD block
eqtl.stat.sig <- eqtl.stat %>% filter(ensg %in% fdr10.ensg$ensg | ensg %in% m2t.stat.sig$ensg)
mqtl.stat.sig <- mqtl.stat %>% filter(cg %in% fdr10.cg$cg)


## show T -> GWAS
n.genes <- length(unique(eqtl.stat.sig$hgnc))
ww = 1 + 5 * 1
hh = 1 + 1 * ceiling(n.genes / 5)

plt <- ggplot(eqtl.stat.sig %>% filter(abs(qtl.z) < 5, abs(gwas.z) < 5), aes(x = qtl.z, y = gwas.z)) +
    theme_bw() + xlab('eQTL') + ylab('GWAS') +
    geom_point(alpha = 0.5, color = '#111199', size = .5) +
    facet_wrap(~hgnc, ncol = 5, scale = 'free')

out.file <- out.hdr %&&% '-genes.pdf'
pdf(file = out.file, useDingbats = FALSE, width = ww, height = hh)
print(plt)
dev.off()

log.msg('%s\n\n', out.file)

## show M -> GWAS
n.cpgs <- length(unique(mqtl.stat.sig$cg))
ww = 1 + 5 * 1
hh = 1 + 1 * ceiling(n.cpgs / 5)

plt <- ggplot(mqtl.stat.sig %>% filter(abs(qtl.z) < 5, abs(gwas.z) < 5), aes(x = qtl.z, y = gwas.z)) +
    theme_bw() + xlab('mQTL') + ylab('GWAS') +
    geom_point(alpha = 0.5, color = '#111199', size = .5) +
    facet_wrap(~cg, ncol = 5, scale = 'free')

out.file <- out.hdr %&&% '-cpgs.pdf'
pdf(file = out.file, useDingbats = FALSE, width = ww, height = hh)
print(plt)
dev.off()

log.msg('%s\n\n', out.file)

## show M -> T
.df <- mqtl.stat.sig %>% filter(cg %in% m2t.stat.sig$cg) %>%
    select(qtl.z, ld, cg, rs) %>% rename(mqtl = qtl.z)

.df <- eqtl.stat.sig %>% filter(ensg %in% m2t.stat.sig$ensg) %>%
    select(qtl.z, ld, hgnc, rs) %>% rename(eqtl = qtl.z) %>%
    left_join(.df, by = c('ld', 'rs')) %>%
    na.omit()

ld.blocks <- unique(.df$ld)

for(.ld in ld.blocks) {

    .df.df <- .df %>% filter(ld == .ld)
    n.genes <- .df.df %>% select(hgnc) %>% unique() %>% nrow()
    n.cpgs <- .df.df %>% select(cg) %>% unique() %>% nrow()

    ww = 1 + n.cpgs
    hh = 1 + n.genes
    plt <- ggplot(.df.df, aes(x = mqtl, y = eqtl)) + theme_bw() +
        geom_point(alpha = 0.5, color = '#111199') +
    facet_grid(hgnc ~ cg, scale = 'free') +
        xlab('mQTL') + ylab('eQTL')

    .ld.str <- gsub(pattern = ':', replacement = '-', .ld)
    out.file <- out.hdr %&&% '-' %&&% .ld.str %&&% '.pdf' 
    pdf(file = out.file, useDingbats = FALSE, width = ww, height = hh)
    print(plt)
    dev.off()

    log.msg('%s\n\n', out.file)
}

log.msg('Done\n\n')
