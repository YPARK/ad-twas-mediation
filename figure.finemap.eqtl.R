argv <- commandArgs(trailingOnly = TRUE)

med.stat.file <- argv[1] # e.g., 'finemap/IGAP_rosmap_eqtl_hs-lm_19.mediation.gz'
qtl.stat.file <- argv[2] # e.g., 'finemap/IGAP_rosmap_eqtl_hs-lm_19.qtl.gz'
out.hdr <- argv[3]

library(dplyr)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = FALSE)
source('figure.util.R')
source('util.R')

med.out.file <- out.hdr %&&% '-global.pdf'
qtl.out.file <- out.hdr %&&% '-qtl.pdf'

qtl.col.names <- c('ensg', 'rs', 'qtl.theta', 'qtl.se', 'qtl.lodds',
                   'chr', 'snp.loc', 'a1', 'a2', 'gwas.theta', 'gwas.se')
qtl.stat <- read.table(qtl.stat.file, col.names = qtl.col.names)

med.col.names <- c('ensg', 'theta', 'theta.var', 'lodds',
                   'theta.2', 'theta.se.2', 'lodds.2',
                   'hgnc', 'tss', 'tes', 'strand',
                   'gwas.theta', 'gwas.z', 'chr', 'ld.start', 'ld.end',
                   'n.snps')

med.stat <- read.table(med.stat.file, col.names = med.col.names) %>%
    arrange(tss) %>% filter(n.snps > 1000) %>%
        mutate(loc = (tss+tes)/2/1e6, theta.se = sqrt(theta.var),
               gwas.lnp = pmin(-log10(2 * pnorm(abs(gwas.z), lower.tail = FALSE)), 20))

med.stat.sig <- med.stat %>% filter(lodds > 2)
qtl.stat <- qtl.stat %>% filter(ensg %in% med.stat.sig$ensg)
med.stat.sig.lab <- med.stat.sig %>% filter(ensg %in% qtl.stat$ensg)

p1 <- ggplot() + theme_classic() + xlab('Genomic Location (Mb)') + ylab('Mediation effect') +
    geom_hline(yintercept = 0, color = 'gray80') +
    geom_point(data = med.stat, aes(x = loc, y = theta), size = 1, color = 'gray60') +
    geom_linerange(data = med.stat.sig,
                   aes(x = loc, ymin = theta - 2*theta.se, ymax = theta + 2*theta.se)) +
    geom_point(data = med.stat.sig, aes(x = loc, y = theta), size = 1, color = 'red') +
    geom_point(data = med.stat.sig.lab, aes(x = loc, y = theta),
               color = 'green', pch = 22, size = 3) +
    geom_text_repel(data = med.stat.sig.lab, aes(x = loc, y = theta, label = hgnc),
                    size = 3)

p0 <- ggplot() + theme_classic() + xlab('Genomic Location (Mb)') +
    ylab('Best GWAS (-log P)') +
    geom_point(data = med.stat, aes(x = loc, y = gwas.lnp), size = 1, color = 'gray') +
    geom_point(data = med.stat.sig, aes(x = loc, y = gwas.lnp), size = 1, color = 'red') +
    geom_point(data = med.stat.sig.lab, aes(x = loc, y = gwas.lnp),
               color = 'green', pch = 22, size = 3) +
    geom_text_repel(data = med.stat.sig.lab, aes(x = loc, y = gwas.lnp, label = hgnc),
                    size = 3)

x.len <- range(med.stat$loc)
x.len <- x.len[2] - x.len[1]

pdf(file = med.out.file, width = x.len/40 + 2, height = 6)
grid.vcat(list(p0, p1), heights = c(2, 2))
dev.off()

################################################################
## 2. for each LD block
for(.ld in unique(med.stat.sig$ld.start)) {

    .stat <- med.stat %>% filter(ld.start == .ld) %>%
        mutate(pr = 1/(1+exp(-lodds))) %>%
            mutate(tss = tss / 1e6, tes = tes / 1e6)

    p0 <- 
        ggplot() + theme_classic() + xlab('Genomic Location (Mb)') + ylab('Best GWAS') +
            geom_point(data = .stat, aes(x = loc, y = gwas.lnp), color = 'gray') +
                geom_point(data = .stat %>% filter(lodds > 0),
                           aes(x = loc, y = gwas.lnp), pch = 22, color = 'orange', size = 3)

    .arrow <- arrow(length = unit(0.0075, 'npc'), angle = 45)

    p1 <-
        ggplot() + theme_classic() + xlab('Genomic Location (Mb)') + ylab('Mediation Prob') +
            geom_point(data = .stat %>% filter(lodds < 0), aes(x = loc, y = pr), color = 'gray')

    p1 <- p1 +
        geom_segment(data = .stat %>% filter(lodds > 0, strand == '+'),
                     aes(x = tss - .01, xend = tes + .01, y = pr, yend = pr),
                     arrow = .arrow, size = 1, color = 'orange')

    p1 <- p1 +
        geom_segment(data = .stat %>% filter(lodds > 0, strand == '-'),
                     aes(x = tes - .01, xend = tss + .01, y = pr, yend = pr),
                     arrow = .arrow, size = 1, color = 'orange')

    p1 <- p1 +
        geom_text_repel(data = .stat %>% filter(lodds > 0),
                        aes(x = loc, y = pr, label = hgnc),
                        size = 3)

    .len <- range(.stat$loc)
    .len <- .len[2] - .len[1]
    .ld.out.file <- out.hdr %&&% '-ld-' %&&% .ld %&&% '.pdf'

    pdf(file = .ld.out.file, width = ceiling(.len) + 1, height = 6)
    print(grid.vcat(list(p0, p1), heights = c(2, 2)))
    dev.off()

}

################################################################
## 3. for each gene
get.qtl.plot <- function(.ensg) {

    .hgnc <- med.stat.sig %>% filter(ensg == .ensg) %>% select(hgnc) %>% unique()
    .hgnc <- as.character(.hgnc[1])

    .stat <- qtl.stat %>% filter(ensg == .ensg) %>%
        select(rs, qtl.theta, qtl.se, gwas.theta, gwas.se)

    .stat.0 <- data.frame(rs = '', qtl.theta = 0, qtl.se = 1, gwas.theta = 0, gwas.se = 1)

    if(nrow(.stat) == 1) {
        .stat <- rbind(.stat, .stat.0)
    }

    ret <- ggplot(.stat, aes(x = qtl.theta/qtl.se, y = gwas.theta/gwas.se, label = rs)) +
        geom_smooth(method = 'lm', se = FALSE) + geom_text_repel(size = 2) +
            geom_point() + xlab('Causal eQTL effect') + ylab('GWAS effect') +
                ggtitle(.hgnc)
    return(ret)
}

ensg.vec <- unique(qtl.stat$ensg)
qtl.plot.list <- lapply(lapply(ensg.vec, get.qtl.plot), ggplotGrob)

ncol <- 5
nrow <- ceiling(length(qtl.plot.list) / ncol)

pdf(file = qtl.out.file, width = 1 + ncol * 1.5, height = nrow + 2)
grid.arrange(grobs = qtl.plot.list, nrow = nrow, ncol = ncol, newpage = TRUE)
dev.off()
