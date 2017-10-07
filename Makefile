CHR := $(shell seq 1 22)
TEMP := /broad/hptmp/ypp/AD/mediation/
QTL_DATA := hs-lm
NULL_DATA := raw-y0

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
# pregenerate mediation data
step2: jobs/step2-jobs.txt.gz
step2-null: jobs/step2-null-jobs.txt.gz

step2-post:   $(foreach chr, $(CHR), \
  stat/IGAP/ld/$(chr).ld.gz \
  $(foreach qtl_data, $(QTL_DATA) $(NULL_DATA), \
  stat/IGAP/data/$(qtl_data)/$(chr).eqtl_bed.gz \
  stat/IGAP/data/$(qtl_data)/$(chr).mqtl_bed.gz ))

jobs/step2-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/obs-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

jobs/step2-null-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/null-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

# % = $(chr)
jobs/step2/obs-%.jobs: gene.batches/%.txt.gz cpg.batches/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA) ; do zcat gene.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP -vQD=$${qtl_data} '{ batch=$$2; print "./make.eqtl.assigned.R" OFS ("eqtl/" CHR "/" batch ".qtl-" QD ".gz") OFS ("eqtl/" CHR "/" batch ".snps.gz") OFS ("eqtl/" CHR "/" batch ".genes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/matched/eqtl/" GWAS "/" QD "/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA) ; do zcat cpg.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP -vQD=$${qtl_data} '{ batch=$$2; print "./make.mqtl.assigned.R" OFS ("mqtl/" CHR "/" batch ".qtl-" QD ".gz") OFS ("mqtl/" CHR "/" batch ".snps.gz") OFS ("mqtl/" CHR "/" batch ".probes.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/matched/mqtl/" GWAS "/" QD "/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

jobs/step2/null-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	zcat gene.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP '{ batch=$$2; print "./make.null.qtl.assigned.R" OFS ("eqtl/" CHR "/" batch ".qtl-$(NULL_DATA).gz") OFS ("eqtl/" CHR "/" batch ".snps.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/null/eqtl/" GWAS "/$(NULL_DATA)/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@
	zcat cpg.batches/$*.txt.gz | awk -F'/' -vOFS=" " -vCHR=$* -vGWAS=IGAP '{ batch=$$2; print "./make.mqtl.assigned.R" OFS ("mqtl/" CHR "/" batch ".qtl-$(NULL_DATA).gz") OFS ("mqtl/" CHR "/" batch ".snps.gz") OFS ("data/" GWAS "/chr" CHR ".ld_bed.gz") OFS ("$(TEMP)/null/mqtl/" GWAS "/$(NULL_DATA)/" CHR "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@

gene.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 eqtl/$*/b*.genes.gz | awk -F'/' '{ gsub(/.genes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | gzip > $@

cpg.batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 mqtl/$*/b*.probes.gz | awk -F'/' '{ gsub(/.probes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | gzip > $@


# % = $(qtl_data)/$(chr)
stat/IGAP/data/%.eqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/eqtl/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@

stat/IGAP/data/%.mqtl_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/mqtl/IGAP/$*/b*.data.gz | xargs zcat) | sort -k1,1 -k2,2g | gzip > $@

# % = $(chr)
stat/IGAP/ld/%.ld.gz : stat/IGAP/data/hs-lm/%.eqtl_bed.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk -F'\t' '{ ld[$$12] ++ } END { for(l in ld) print l FS ld[l] }' | sed 's/:/\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > $@


################################################################
# Fine-mapping and bootstrapping
step3: jobs/step3-eqtl-jobs.txt.gz \
  jobs/step3-mqtl-jobs.txt.gz #\
  jobs/step3-m2t-jobs.txt.gz

step3-resubmit: jobs/step3-eqtl-jobs-resubmit.txt.gz jobs/step3-mqtl-jobs-resubmit.txt.gz # jobs/step3-m2t-jobs-resubmit.txt.gz

jobs/step3-%-jobs-resubmit.txt.gz: jobs/step3-%-jobs.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0 && system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Re-$*-ZQTL -binding "linear:1" -l h_rt=7200 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/step3-%-jobs.txt.gz: $(foreach chr, $(CHR), $(foreach task, nwas, jobs/step3/$(task)-%-$(chr).jobs)) # finemap bootstrap 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@cat $^ | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N $*-ZQTL -binding "linear:1" -l h_rt=1000 -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

jobs/step3/finemap-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-finemap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS ("finemap/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done

jobs/step3/bootstrap-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-bootstrap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS "direct" FS ("bootstrap/direct/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-bootstrap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS "marginal" FS ("bootstrap/marginal/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done

jobs/step3/nwas-eqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.nwas.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS "TRUE" FS ("nwas/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR ".nwas.gz") }' >> $@; done

jobs/step3/nwas-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.nwas.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS "FALSE" FS ("nwas/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR ".nwas.gz") }' >> $@; done

jobs/step3/finemap-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-finemap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS ("finemap/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done

jobs/step3/bootstrap-mqtl-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-bootstrap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS "direct" FS ("bootstrap/direct/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-bootstrap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS "marginal" FS ("bootstrap/marginal/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done

jobs/step3/finemap-m2t-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-m2t.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS 356 OFS ("m2t/" GWAS "/rosmap/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done


step3-post: $(foreach data, $(QTL_DATA), $(foreach chr, $(CHR), $(foreach stat, mediation qtl, m2t/IGAP_rosmap_$(data)_$(chr).$(stat).gz \
 $(foreach reg, eqtl mqtl, finemap/IGAP_rosmap_$(reg)_$(data)_$(chr).$(stat).gz) ) ) ) \
 $(foreach data, $(QTL_DATA), $(foreach null, direct marginal, $(foreach reg, eqtl mqtl, $(foreach chr, $(CHR), bootstrap/$(null)_IGAP_rosmap_$(reg)_hs-lm_$(chr).mediation.gz ) ) ) ) \
 $(foreach data, $(QTL_DATA), $(foreach reg, eqtl mqtl, $(foreach chr, $(CHR), nwas/IGAP_rosmap_$(reg)_hs-lm_$(chr).nwas.gz ) ) )

# % = IGAP_rosmap_eqtl_hs-lm_21
finemap/%.mediation.gz:
	(ls -1 finemap/$(shell echo $* | sed 's/_/\//g')/*.mediation.gz 2> /dev/null | xargs cat) > $@

finemap/%.qtl.gz:
	(ls -1 finemap/$(shell echo $* | sed 's/_/\//g')/*.qtl.gz 2> /dev/null | xargs zcat) | gzip > $@

m2t/%.mediation.gz:
	(ls -1 m2t/$(shell echo $* | sed 's/_/\//g')/*.mediation.gz 2> /dev/null | xargs cat) > $@

m2t/%.qtl.gz:
	(ls -1 m2t/$(shell echo $* | sed 's/_/\//g')/*.qtl.gz 2> /dev/null | xargs cat) > $@

# % = $(null)_IGAP_rosmap_eqtl_hs-lm_$(chr)
bootstrap/%.mediation.gz:
	(ls -1 bootstrap/$(shell echo $* | sed 's/_/\//g')/*.mediation.gz 2> /dev/null | xargs zcat) | awk 'NF > 0' | gzip > $@

nwas/%.nwas.gz:
	(ls -1 nwas/$(shell echo $* | sed 's/_/\//g')/*.nwas.gz 2> /dev/null | xargs zcat) | awk 'NF > 0' | gzip > $@

################################################################
step3-figure: jobs/step3-bootstrap-figures.jobs.gz

jobs/step3-bootstrap-figures.jobs.gz: tables/bootstrap_ld_strict.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	awk -vN=$$(zcat $< | tail -n+2 | wc -l) 'BEGIN { for(j=1; j<=N; ++j) print "./figure.bootstrap.ld.R" FS j FS  "figures/ld/bootstrap/" }' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N FIG-ZQTL -binding "linear:1" -l h_rt=5000 -l h_vmem=16g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

## finemap/figures/IGAP_rosmap_%-global.pdf: finemap/IGAP_rosmap_eqtl_hs-lm_%.mediation.gz
##	mkdir -p $(dir $@)
##	Rscript --vanilla figure.finemap.R finemap/IGAP_rosmap_eqtl_hs-lm_$*.mediation.gz finemap/IGAP_rosmap_mqtl_hs-lm_$*.mediation.gz m2t/IGAP_rosmap_hs-lm_$*.mediation.gz finemap/figures/IGAP_rosmap_$*

step3-table: tables/bootstrap_gene_significant.txt.gz

tables/bootstrap_ld_strict.txt.gz: table.bootstrap.R
	Rscript --vanilla $<


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

