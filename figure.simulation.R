
options(stringsAsFators = FALSE)
source('util.R')
source('figure.util.R')
library(dplyr)
library(PRROC)
library(broom)
library(ggplot2)
library(readr)

get.auprc <- function(score, lab) {
    pr.out <- pr.curve(score[lab == 1], score[lab == 0])
    return(pr.out$auc.davis.goadrich)
}

get.power <- function(score, lab, fdr.cutoff = 0.01) {
    pr.out <- pr.curve(score[lab == 1], score[lab == 0], curve = TRUE)
    .idx <- pr.out$curve[, 2] >= (1 - fdr.cutoff)
    if(sum(.idx) == 0) return(0)
    .temp <- pr.out$curve[.idx, 1]
    return(.temp[1])
}

zqtl.power <- function(tab) {
    as.numeric(get.power(tab$lodds, tab$causal, 0.01))
}

twas.power <- function(tab) {
    as.numeric(get.power(abs(tab$twas), tab$causal, 0.01))
}

egger.power <- function(tab) {
    as.numeric(get.power(abs(tab$egger.t), tab$causal, 0.01))
}

zqtl.auprc <- function(tab) {
    as.numeric(get.auprc(tab$lodds, tab$causal))
}

twas.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$twas), tab$causal))
}

egger.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$egger.t), tab$causal))
}

summarize.sim.tab <- function(sim.tab, egger.score, twas.score, zqtl.score) {

    temp.egger <- sim.tab %>%
        group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med, n.med) %>%
            do(egger = egger.score(.)) %>% tidy(egger) %>%
                group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med) %>%
                    summarize(score = mean(x), score.se = sd(x)/sqrt(length(x)))
    
    temp.twas <- sim.tab %>%
        group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med, n.med) %>%
            do(twas = twas.score(.)) %>% tidy(twas) %>%
                group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med) %>%
                    summarize(score = mean(x), score.se = sd(x)/sqrt(length(x)))
    
    temp.zqtl <- sim.tab %>%
        group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med, n.med) %>%
            do(zqtl = zqtl.score(.)) %>% tidy(zqtl) %>%
                group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med) %>%
                    summarize(score = mean(x), score.se = sd(x)/sqrt(length(x)))
    
    .df <- rbind(temp.zqtl %>% mutate(method = 'zqtl'),
                 temp.twas %>% mutate(method = 'twas'),
                 temp.egger %>% mutate(method = 'egger'))
    
    .df$method <- factor(.df$method, c('zqtl', 'twas', 'egger'))
    .df$n.causal.med <- factor(.df$n.causal.med, 1:3, 'causal mediators = ' %&&% 1:3)
    .df$pve.dir <- factor(.df$pve.dir, unique(.df$pve.dir), 'pleiotropy = ' %&&% unique(.df$pve.dir))

    return(.df)
}

plt.sim <- function(.df, n.qtl) {

    plt.df <- .df %>% filter(n.causal.qtl == n.qtl)
    plt.aes <- aes(x=pve.med, y = score, color = method, shape = method,
                   ymin = score - 2 * score.se, ymax = score + 2 * score.se)

    plt <- ggplot(plt.df, plt.aes) + theme_bw() + xlab('Mediation PVE') +
        geom_errorbar(width = .005) + geom_point(size = 2, fill = 'white') +
            geom_line() +
                facet_grid(n.causal.med ~ pve.dir, scale = 'free') +
                    ggtitle('QTLs per gene = ' %&&% n.qtl) +
                        theme(legend.title = element_blank(),
                              axis.text.x = element_text(angle=45))
    
    plt <- plt + scale_shape_manual(values = c(21, 22, 24))

    return(plt)
}

################################################################
## Simulation 1

ld.sim.cols <- c('gene', 'causal', 'pve.med', 'pve.dir', 'pve.qtl',
                 'p', 'n.med', 'n.causal.qtl', 'n.causal.med', 'n.causal.direct',
                 'theta', 'theta.se', 'lodds',
                 'egger.t', 'twas')

sim.files <- 'simulation/result/ld_IGAP_rosmap_eqtl_hs-lm_' %&&% (1:22) %&&% '.sim.gz'

sim.tab <- do.call(rbind, lapply(sim.files, read_tsv, col_names = ld.sim.cols))

sim.tab.power <- summarize.sim.tab(sim.tab, egger.power, twas.power, zqtl.power)

sim.tab.auprc <- summarize.sim.tab(sim.tab, egger.auprc, twas.auprc, zqtl.auprc)

for(q in 1:3) {

    pdf(file = 'figures/simulation1_power_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = FALSE)
    print(plt.sim(sim.tab.power, n.qtl = 3) + ylab('Power (FDR < 1%)'))
    dev.off()

    pdf(file = 'figures/simulation1_auprc_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = FALSE)
    print(plt.sim(sim.tab.auprc, n.qtl = 3) + ylab('AUPRC'))
    dev.off()

}

################################################################
## Simulation 2

pleio.sim.cols <- c('gene', 'causal', 'pleiotropy', 'pve.med', 'pve.dir', 'pve.qtl',
                    'p', 'n.med', 'n.causal.qtl', 'n.causal.med', 'n.causal.direct',
                    'theta', 'theta.se', 'lodds',
                    'egger.t', 'twas')

sim.files <- 'simulation/result/pleiotropy_IGAP_rosmap_eqtl_hs-lm_' %&&% (1:22) %&&% '.sim.gz'
sim.tab <- do.call(rbind, lapply(sim.files, read_tsv, col_names = pleio.sim.cols))

zqtl.false.power <- function(tab) {
    as.numeric(get.power(tab$lodds, -tab$pleiotropy, 0.5))
}

twas.false.power <- function(tab) {
    as.numeric(get.power(abs(tab$twas), -tab$pleiotropy, 0.5))
}

egger.false.power <- function(tab) {
    as.numeric(get.power(abs(tab$egger.t), -tab$pleiotropy, 0.5))
}

zqtl.false.auprc <- function(tab) {
    as.numeric(get.auprc(tab$lodds, -tab$pleiotropy))
}

twas.false.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$twas), -tab$pleiotropy))
}

egger.false.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$egger.t), -tab$pleiotropy))
}

################################################################
## negative prediction
sim.false.tab.power <-
    summarize.sim.tab(sim.tab, egger.false.power, twas.false.power, zqtl.false.power)

sim.false.tab.auprc <-
    summarize.sim.tab(sim.tab, egger.false.auprc, twas.false.auprc, zqtl.false.auprc)

for(q in 1:3) {

    pdf(file = 'figures/simulation2_negative_power_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = FALSE)
    print(plt.sim(sim.false.tab.power, n.qtl = 3) + ylab('Mis-classification Power (FDR < 50%)'))
    dev.off()

    pdf(file = 'figures/simulation2_negative_auprc_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = FALSE)
    print(plt.sim(sim.false.tab.auprc, n.qtl = 3) + ylab('Mis-classification AUPRC'))
    dev.off()

}

################################################################
## positive prediction
sim.true.tab.power <-
    summarize.sim.tab(sim.tab, egger.power, twas.power, zqtl.power)

sim.true.tab.auprc <-
    summarize.sim.tab(sim.tab, egger.auprc, twas.auprc, zqtl.auprc)

for(q in 1:3) {
    pdf(file = 'figures/simulation2_positive_power_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = TRUE)
    print(plt.sim(sim.true.tab.power, n.qtl = 3) + ylab('(FDR < 1%)'))
    dev.off()

    pdf(file = 'figures/simulation2_positive_auprc_nqtl' %&&% q %&&% '.pdf', width = 8, height = 5, useDingbats = TRUE)
    print(plt.sim(sim.true.tab.auprc, n.qtl = 3) + ylab('AUPRC'))              
    dev.off()
}
