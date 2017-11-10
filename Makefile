CHR := $(shell seq 1 22)
TEMP := /broad/hptmp/ypp/AD/mediation/
QTL_DATA := hs-lm
NULL_DATA := raw-y0
TIS := CO HC

all:

################################################################
# combine LD blocks
step1: $(foreach chr, $(CHR), data/IGAP/chr$(chr).ld_bed.gz)

# % = $(chr)
data/IGAP/chr%.ld_bed.gz: LD/igap-ad.bed.gz IGAP/chr%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@./make_ld_blocks.sh $* $^ $@

$(TEMP) :
	[ -d $@ ] || mkdir -p $@

################################################################
# Pre-generate mediation data
QTL_CUTOFF := 2

step2: jobs/step2-jobs.txt.gz # jobs/step2-hic-jobs.txt.gz

step2-long: jobs/step2-jobs-long.txt.gz # jobs/step2-hic-jobs-long.txt.gz

step2-post-qtl := $(foreach chr, $(CHR), stat/IGAP/ld/$(chr).ld.gz \
  $(foreach qtl_data, $(QTL_DATA), \
  stat/IGAP/data/$(qtl_data)/$(chr).mqtl_bed.gz \
  stat/IGAP/data/$(qtl_data)/$(chr).eqtl_bed.gz ))

step2-post-hic := $(foreach chr, $(CHR), \
  $(foreach qd, $(QTL_DATA), \
  stat/IGAP/hic-merged/$(qd)/$(chr).eqtl_bed.gz \
  stat/IGAP/hic-merged/$(qd)/$(chr).mqtl_bed.gz \
  $(foreach tis, $(TIS), \
  stat/IGAP/hic-data/$(tis)/$(qd)/$(chr).eqtl_bed.gz \
  stat/IGAP/hic-data/$(tis)/$(qd)/$(chr).mqtl_bed.gz )))

step2-post: $(step2-post-qtl) $(step2-post-hic)

## step2-queue:
##	@[ $$(zcat jobs/step2-jobs.txt.gz | wc -l) -lt 1 ] || qsub -t 1-$$(zcat jobs/step2-jobs.txt.gz | wc -l) -N med.data -binding "linear:1" -l h_rt=1800 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh jobs/step2-jobs.txt.gz

jobs/step2-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/obs-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.data -binding "linear:1" -l h_rt=1800 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

# % = $(chr)
jobs/step2/obs-%.jobs: gene.batches/%.txt.gz cpg.batches/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA) ; do zcat gene.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP -vQD=$${qtl_data} '{ batch=$$2; print "./make.eqtl.assigned.R" OFS ("eqtl/" CHR "/" batch ".qtl-" QD ".gz") OFS ("eqtl/" CHR "/" batch ".snps.gz") OFS ("eqtl/" CHR "/" batch ".genes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/matched/eqtl/" GWAS "/" QD "/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA) ; do zcat cpg.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP -vQD=$${qtl_data} '{ batch=$$2; print "./make.mqtl.assigned.R" OFS ("mqtl/" CHR "/" batch ".qtl-" QD ".gz") OFS ("mqtl/" CHR "/" batch ".snps.gz") OFS ("mqtl/" CHR "/" batch ".probes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/matched/mqtl/" GWAS "/" QD "/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

gene.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 eqtl/$*/b*.genes.gz | awk -F'/' '{ gsub(/.genes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | gzip > $@

cpg.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 mqtl/$*/b*.probes.gz | awk -F'/' '{ gsub(/.probes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | gzip > $@


# % = $(qtl_data)/$(chr)
stat/IGAP/data/%.eqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/eqtl/IGAP/$*/b*.data.gz | xargs zcat) | awk '$$8 > $(QTL_CUTOFF) || $$8 < -$(QTL_CUTOFF)' | sort -k1,1 -k2,2g | gzip > $@

stat/IGAP/data/%.mqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/mqtl/IGAP/$*/b*.data.gz | xargs zcat) | awk '$$8 > $(QTL_CUTOFF) || $$8 < -$(QTL_CUTOFF)' | sort -k1,1 -k2,2g | gzip > $@

# % = $(chr)
stat/IGAP/ld/%.ld.gz : stat/IGAP/data/hs-lm/%.eqtl_bed.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk -F'\t' '{ ld[$$12] ++ } END { for(l in ld) print l FS ld[l] }' | sed 's/:/\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > $@

jobs/step2-jobs-long.txt.gz: jobs/step2-jobs.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.hic.data-long -binding "linear:1" -l h_rt=7200 -l h_vmem=8g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

################################################################
# Hi-C guided filtering of QTL pairs
jobs/step2-hic-jobs-long.txt.gz: jobs/step2-hic-jobs.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.hic.data-long -binding "linear:1" -l h_rt=7200 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step2-hic-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/hic-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.hic.data -binding "linear:1" -l h_rt=1800 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step2/hic-%.jobs: gene-hic.batches/%.txt.gz cpg-hic.batches/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@

	for tis in $(TIS); do for qtl_data in $(QTL_DATA) ; do zcat cpg-hic.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vTIS=$${tis} -vCHR=$* -vGWAS=IGAP -vQT=mqtl -vQD=$${qtl_data} '{ batch=$$2; IN = ("hic-" QT "/" TIS "/" CHR "/" batch); OUT = ("$(TEMP)/matched/hic-" QT "/" TIS "/" GWAS "/" QD "/" CHR); print ("./make." QT ".assigned.R") OFS (IN ".qtl-" QD ".gz") OFS (IN ".snps.gz") OFS (IN ".probes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS (OUT "/" batch ".data.gz") }' >> $@; done; done

	for tis in $(TIS); do for qtl_data in $(QTL_DATA) ; do zcat gene-hic.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vTIS=$${tis} -vCHR=$* -vGWAS=IGAP -vQT=eqtl -vQD=$${qtl_data} '{ batch=$$2; IN = ("hic-" QT "/" TIS "/" CHR "/" batch); OUT = ("$(TEMP)/matched/hic-" QT "/" TIS "/" GWAS "/" QD "/" CHR); print ("./make." QT ".assigned.R") OFS (IN ".qtl-" QD ".gz") OFS (IN ".snps.gz") OFS (IN ".genes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS (OUT "/" batch ".data.gz") }' >> $@; done; done

gene-hic.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 hic-eqtl/*/$*/b*.genes.gz | awk -F'/' '{ gsub(/.genes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | awk '{ uniq[$$0]++ } END { for(u in uniq) print u }' | sort | gzip > $@

cpg-hic.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 hic-mqtl/*/$*/b*.probes.gz | awk -F'/' '{ gsub(/.probes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | awk '{ uniq[$$0]++ } END { for(u in uniq) print u }' | sort | gzip > $@


################################################################
# % = $(qtl_data)/$(chr)
stat/IGAP/hic-merged/%.eqtl_bed.gz: $(foreach tis, $(TIS), stat/IGAP/hic-data/$(tis)/%.eqtl_bed.gz)
	./run.sh ./make.merge.qtl.R TRUE $@ $^

stat/IGAP/hic-merged/%.mqtl_bed.gz: $(foreach tis, $(TIS), stat/IGAP/hic-data/$(tis)/%.mqtl_bed.gz)
	./run.sh ./make.merge.qtl.R FALSE $@ $^

stat/IGAP/hic-data/CO/%.eqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/hic-eqtl/CO/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@

stat/IGAP/hic-data/CO/%.mqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/hic-mqtl/CO/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@

stat/IGAP/hic-data/HC/%.eqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/hic-eqtl/HC/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@

stat/IGAP/hic-data/HC/%.mqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/hic-mqtl/HC/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@


################################################################
# Model estimation, bootstrap, PVE, etc.
step3: jobs/step3-eqtl-jobs.txt.gz jobs/step3-mqtl-jobs.txt.gz

step3-resubmit: jobs/step3-eqtl-jobs-resubmit.txt.gz jobs/step3-mqtl-jobs-resubmit.txt.gz

jobs/step3-%-jobs-resubmit.txt.gz: jobs/step3-%-jobs.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0 && system("[ ! -f " $$NF " ]") == 0 && system("[ ! -f " $$NF ".gene-mediation.gz ]") == 0 && system("[ ! -f " $$NF ".cpg-mediation.gz ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Re-$*-ZQTL -binding "linear:1" -l h_rt=24:00:00 -l h_vmem=8g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step3-%-jobs.txt.gz: $(foreach chr, $(CHR), $(foreach task, full hic, jobs/step3/$(task)-%-$(chr).jobs))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@cat $^ | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N $*-ZQTL -binding "linear:1" -l h_rt=1800 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step3/nwas-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.nwas.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS "TRUE" FS ("nwas/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR ".nwas.gz") }' >> $@; done

jobs/step3/nwas-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.nwas.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS "FALSE" FS ("nwas/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR ".nwas.gz") }' >> $@; done

jobs/step3/hic-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/hic-merged/" QD "/" chr ".eqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") OFS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS Z FS ("result/hic-" Z "/" GWAS "/eqtl/" QD "/" chr "/" NR) }' >> $@; done; done

jobs/step3/hic-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/hic-merged/" QD "/" chr ".mqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS Z FS ("result/hic-" Z "/" GWAS "/mqtl/" QD "/" chr "/" NR) }' >> $@; done; done

jobs/step3/hic-joint-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-joint.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/hic-merged/" QD "/" chr ".mqtl_bed.gz") OFS ("stat/" GWAS "/hic-merged/" QD "/" chr ".eqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS Z FS ("result/hic-" Z "/" GWAS "/joint/" QD "/" chr "/" NR) }' >> $@; done; done

jobs/step3/full-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS Z FS ("result/full-" Z "/" GWAS "/eqtl/" QD "/" chr "/" NR) }' >> $@; done; done

jobs/step3/full-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS Z FS ("result/full-" Z "/" GWAS "/mqtl/" QD "/" chr "/" NR) }' >> $@; done; done

jobs/step3/full-joint-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for z_cutoff in 5 3 0; do for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vZ=$${z_cutoff} -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-joint.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") OFS ("IGAP/chr" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS Z FS ("result/full-" Z "/" GWAS "/joint/" QD "/" chr "/" NR) }' >> $@; done; done

################################################################
# combine all data
step3-post : $(foreach inter, full hic, $(foreach cutoff, 5 3 0, $(foreach qtl, eqtl, $(foreach chr, $(CHR), $(foreach stat, mediation null, result/$(stat)/$(inter)-$(cutoff)_IGAP_$(qtl)_hs-lm_$(chr).$(stat).gz)))))

# % = $(hic-8)_IGAP_$(mqtl)_hs-lm_$(21)
result/mediation/%.mediation.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 result/$(shell echo $* | sed 's/_/\//g')/*.mediation.gz 2> /dev/null | xargs zcat) | awk 'NF > 0' | gzip > $@

# % = $(hic-8)_IGAP_$(mqtl)_hs-lm_$(21)
result/null/%.null.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 result/$(shell echo $* | sed 's/_/\//g')/*.null.gz 2> /dev/null | xargs zcat) | awk 'NF > 0' | gzip > $@


################################################################
step3-figure: jobs/step3-figures.jobs.gz

jobs/step3-figures.jobs.gz: $(foreach d, $(shell ls -1 tables/genes/ 2> /dev/null), jobs/step3/$(d)-figures.jobs.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N figure -binding "linear:1" -l h_rt=3600 -l h_vmem=8g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step3/%-figures.jobs.gz: tables/genes/%/significant_LD.txt.gz 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" | gzip > $@
	zcat $< | tail -n+2 | awk '{ k=$$1 FS $$2 FS $$3; pve = $$(NF - 4); if(!(k in keys) || keys[k] < pve) { keys[k] = pve; gene[k] = $$8} } END { for(k in keys) print k FS keys[k] FS gene[k] }'| sort -k1 -k2n | awk '{ printf "./make.gene.figure.local.R %s %d %d %d figures/genes/local_$*/local_chr%d_%0.0f_%0.0f_%s.pdf\n", "$<", $$1, $$2, $$3, $$1, ($$2/1e6), ($$3/1e6), $$5 }' | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip >> $@


################################################################
## simulation analysis
step4: $(foreach chr, $(CHR), simulation/large-ld-$(chr).txt.gz) \
       jobs/step4-eqtl-jobs.txt.gz

jobs/step4-%-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step4/simulation-%-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@cat $^ | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N $*-SIM -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

step4-resubmit: jobs/step4-eqtl-jobs-resubmit.txt.gz

jobs/step4-%-jobs-resubmit.txt.gz: jobs/step4-%-jobs.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Re-$*-SIM -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

simulation/large-ld-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh ./make.simulation-select.R stat/IGAP/ld/$*.ld.gz geno/rosmap1709-chr$* $@

jobs/step4/simulation-eqtl-%.jobs: simulation/large-ld-%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat simulation/large-ld-$*.txt.gz | awk '{ print $$0 FS NR }' | shuf | head -n3 | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} 'BEGIN { jobid = 0 } { prog = "./make.simulation.R" FS ("simulation/large-ld-" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS $$NF; rdir = "simulation/ld/" GWAS "/rosmap/eqtl/" QD "/" chr; m = split("0.01,0.02,0.05,0.1,0.15,0.2",h1_arr,","); for(j=1; j<=m; ++j) { h1=h1_arr[j]; for(h2=0; h2<=.2; h2+=.05) { for(c1=1; c1<=3; c1++) { for(c2=1; c2<=3; ++c2) { print prog FS 0.17 FS h1 FS h2 FS c1 FS c2 FS c2 FS (rdir "/sim-" (++jobid) ".gz") } } } } }' >> $@; done
	for qtl_data in $(QTL_DATA); do zcat simulation/large-ld-$*.txt.gz | awk '{ print $$0 FS NR }' | shuf | head -n3 | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} 'BEGIN { jobid = 0 } { prog = "./make.simulation-pleiotropy.R" FS ("simulation/large-ld-" chr ".txt.gz") FS ("geno/rosmap1709-chr" chr) FS $$NF; rdir = "simulation/pleiotropy/" GWAS "/rosmap/eqtl/" QD "/" chr; m = split("0.01,0.02,0.05,0.1,0.15,0.2",h1_arr,","); for(j=1; j<=m; ++j) { h1=h1_arr[j]; for(h2=0; h2<=.2; h2+=.05) { for(c1=1; c1<=3; c1++) { for(c2=1; c2<=3; ++c2) { print prog FS 0.17 FS h1 FS h2 FS c1 FS c2 FS c2 FS (rdir "/sim-" (++jobid) ".gz") } } } } }' >> $@; done

step4-post: $(foreach sim, pleiotropy ld, $(foreach chr, $(shell seq 1 22), simulation/result/$(sim)_IGAP_rosmap_eqtl_hs-lm_$(chr).sim.gz))

simulation/result/%.sim.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 simulation/$(shell echo $* | sed 's/_/\//g')/sim-*.gz 2> /dev/null | xargs zcat) | gzip > $@


################################################################
## Utilities
PLINKZIP := https://www.cog-genomics.org/static/bin/plink170906/plink_linux_x86_64.zip

bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/
