#!/usr/bin/env Rscript

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

get.pleiotropy.f1 <- function(score, lab, lab.false) {

    ## The pr curve as a matrix, where the first column contains recall,
    ## the second contains precision, and the third contains the
    ## corresponding threshold on the scores.
    ##
    ## Take argmax F1 score threshold and determine FDR due to pleiotropy
    ##

    pr.out <- pr.curve(score[lab == 1], score[lab == 0], curve = TRUE)

    pr.curve <- pr.out$curve %>% as.data.frame()
    colnames(pr.curve) <- c('recall', 'precision', 'threshold')

    best.pr <- pr.curve %>%
        mutate(f1 = 2 * (recall * precision)/(recall + precision)) %>%
            slice(which.max(f1))

    ret <- mean(lab.false[score >= best.pr$threshold])
    ## denom <- sum(lab[score >= best.pr$threshold])
    ## ret <- num / pmax(denom, 1)
    return(ret)
}

get.pleiotropy.recall <- function(score, lab, lab.false, recall.cutoff = .2) {

    ## The pr curve as a matrix, where the first column contains recall,
    ## the second contains precision, and the third contains the
    ## corresponding threshold on the scores.
    ##
    ## Take argmax F1 score threshold and determine FDR due to pleiotropy
    ##

    pr.out <- pr.curve(score[lab == 1], score[lab == 0], curve = TRUE)

    pr.curve <- pr.out$curve %>% as.data.frame()
    colnames(pr.curve) <- c('recall', 'precision', 'threshold')

    best.pr <- pr.curve %>% filter(recall >= recall.cutoff) %>%
        slice(which.max(precision))

    ret <- mean(lab.false[score >= best.pr$threshold])
    ## denom <- sum(lab[score >= best.pr$threshold])
    ## ret <- num / pmax(denom, 1)
    return(ret)
}

cammel.power <- function(tab) {
    as.numeric(get.power(tab$lodds, tab$causal, 0.01))
}

twas.power <- function(tab) {
    as.numeric(get.power(abs(tab$twas), tab$causal, 0.01))
}

egger.power <- function(tab) {
    as.numeric(get.power(abs(tab$egger.t), tab$causal, 0.01))
}

cammel.auprc <- function(tab) {
    as.numeric(get.auprc(tab$lodds, tab$causal))
}

twas.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$twas), tab$causal))
}

egger.auprc <- function(tab) {
    as.numeric(get.auprc(abs(tab$egger.t), tab$causal))
}

cammel.false.f1 <- function(tab) {
    as.numeric(get.pleiotropy.f1(tab$lodds, tab$causal, -tab$pleiotropy))
}

twas.false.f1 <- function(tab) {
    as.numeric(get.pleiotropy.f1(abs(tab$twas), tab$causal, -tab$pleiotropy))
}

egger.false.f1 <- function(tab) {
    as.numeric(get.pleiotropy.f1(abs(tab$egger.t), tab$causal, -tab$pleiotropy))
}

cammel.false.25recall <- function(tab) {
    as.numeric(get.pleiotropy.recall(tab$lodds, tab$causal, -tab$pleiotropy, 0.25))
}

twas.false.25recall <- function(tab) {
    as.numeric(get.pleiotropy.recall(abs(tab$twas), tab$causal, -tab$pleiotropy, 0.25))
}

egger.false.25recall <- function(tab) {
    as.numeric(get.pleiotropy.recall(abs(tab$egger.t), tab$causal, -tab$pleiotropy, 0.25))
}

summarize.sim.tab <- function(sim.tab, egger.score, twas.score, cammel.score) {

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

    temp.cammel <- sim.tab %>%
        group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med, n.med) %>%
            do(cammel = cammel.score(.)) %>% tidy(cammel) %>%
                group_by(pve.med, pve.dir, pve.qtl, n.causal.qtl, n.causal.med) %>%
                    summarize(score = mean(x), score.se = sd(x)/sqrt(length(x)))

    .df <- rbind(temp.cammel %>% mutate(method = 'cammel'),
                 temp.twas %>% mutate(method = 'twas'),
                 temp.egger %>% mutate(method = 'egger'))

    .df$method <- factor(.df$method, c('cammel', 'twas', 'egger'))
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

pleio.sim.cols <- c('gene', 'causal', 'pleiotropy', 'pve.med', 'pve.dir', 'pve.qtl',
                    'p', 'n.med', 'n.causal.qtl', 'n.causal.med', 'n.causal.direct',
                    'theta', 'theta.se', 'lodds',
                    'egger.t', 'twas')

sim.files <- 'simulation/result/pleiotropy_IGAP_rosmap_eqtl_hs-lm_' %&&% (1:22) %&&% '.sim.gz'
sim.tab <- do.call(rbind, lapply(sim.files, read_tsv, col_names = pleio.sim.cols))

for(nm in 1:3) {
    for(nq in 1:3) {

        sim.auprc.tab <- 
            summarize.sim.tab(sim.tab %>% filter(pve.dir > 0, n.causal.med == nm, n.causal.qtl == nq),
                              egger.auprc, twas.auprc, cammel.auprc)

        sim.power.tab <- 
            summarize.sim.tab(sim.tab %>% filter(pve.dir > 0, n.causal.med == nm, n.causal.qtl == nq),
                              egger.power, twas.power, cammel.power)

        sim.false.tab <- 
            summarize.sim.tab(sim.tab %>% filter(pve.dir > 0, n.causal.med == nm, n.causal.qtl == nq),
                              egger.false.f1, twas.false.f1, cammel.false.f1)


        .aes <- aes(x = pve.med * 100, y = score, color = method, shape = method,
                    ymin = (score - 2*score.se), ymax = (score + 2*score.se))

        p1 <- gg.plot(sim.auprc.tab, .aes) +
            geom_errorbar(width = .005) + geom_point(size = 2, fill = 'white') +
                geom_line() +
                    facet_grid(n.causal.med~pve.dir) +
                        scale_shape_manual(values = c(21, 22, 24)) +
                            ylab('AUPRC') +
                                xlab('% variance explained by mediation') +
                                    ggtitle('# causal QTLs = ' %&&% nq) +
        theme(legend.position = 'none')

        p2 <- gg.plot(sim.power.tab, .aes) +
            geom_errorbar(width = .005) + geom_point(size = 2, fill = 'white') +
                geom_line() +
                    facet_grid(n.causal.med~pve.dir) +
                        scale_shape_manual(values = c(21, 22, 24)) +
                            ylab('power at FDR < 1%') +
                                xlab('% variance explained by mediation') +
                                    theme(legend.position = 'none')

        .aes <- aes(x = pve.med * 100, y = score * 100, color = method, shape = method,
                    ymin = 100*(score - 2*score.se), ymax = 100*(score + 2*score.se))

        p3 <- gg.plot(sim.false.tab, .aes) +
            geom_errorbar(width = .005) + geom_point(size = 2, fill = 'white') +
                geom_line() +
                    facet_grid(n.causal.med~pve.dir) +
                        scale_shape_manual(values = c(21, 22, 24)) +
                            ylab('% false discovery at max F1') +
                                xlab('% variance explained by mediation') +
                                    theme(legend.position = 'bottom')

        gg <- grid.vcat(list(p1, p2, p3), heights = c(1.2, 1, 1.3))

        out.file <- 'figures/Fig_simulation_med' %&&% nm %&&% '_qtl' %&&% nq %&&% '.pdf'

        ggsave(filename = out.file, plot = gg, width = 8, height = 9)
    }
}
