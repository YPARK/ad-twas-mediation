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
# Fine-mapping
step3: jobs/step3-jobs.txt.gz

step3-queue:
	@[ $$(zcat jobs/step3-jobs.txt.gz | wc -l) -lt 1 ] || qsub -t 1-$$(zcat jobs/step3-jobs.txt.gz | wc -l) -N FM-ZQTL -binding "linear:1" -q short -l h_vmem=3g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh jobs/step3-jobs.txt.gz

step3-resubmit:
	zcat jobs/step3-jobs.txt.gz | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' | gzip > jobs/step3-jobs-resubmit.txt.gz
	@[ $$(zcat jobs/step3-jobs-resubmit.txt.gz | wc -l) -lt 1 ] || qsub -t 1-$$(zcat jobs/step3-jobs-resubmit.txt.gz | wc -l) -N FM-ZQTL -binding "linear:1" -q long -l h_vmem=3g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh jobs/step3-jobs-resubmit.txt.gz

jobs/step3-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step3/finemap-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N FM-ZQTL -binding "linear:1" -q short -l h_vmem=3g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

jobs/step3/finemap-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-finemap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 356 FS "TRUE" FS ("finemap/" GWAS "/rosmap/eqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-finemap.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS "FALSE" FS ("finemap/" GWAS "/rosmap/mqtl/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done
	for qtl_data in $(QTL_DATA); do zcat stat/IGAP/ld/$*.ld.gz | awk -vchr=$* -vGWAS=IGAP -vQD=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation-m2t.R" OFS ("stat/" GWAS "/ld/" chr ".ld.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".mqtl_bed.gz") OFS ("stat/" GWAS "/data/" QD "/" chr ".eqtl_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74046 FS 598 FS 356 OFS ("m2t/" GWAS "/rosmap/" QD "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".mediation.gz ]") == 0' >> $@; done




################################################################
## Utilities
PLINKZIP := https://www.cog-genomics.org/static/bin/plink170906/plink_linux_x86_64.zip

bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/

