
CHR := $(shell seq 1 22)
TEMP := /broad/hptmp/ypp/AD/twas/mediation/
QTL_DATA := hs-fqtl hs-sqtl peer-sqtl hs-lm peer-lm

# BATCHES := $(shell ls -1 qtl/*/b*.genes.gz | awk -F'/' '{ gsub(/.genes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }')

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
step2-perm: jobs/step2-perm-jobs.txt.gz

step2-post: $(foreach qtl_data, $(QTL_DATA), \
  $(foreach chr, $(CHR), \
  $(foreach stat, obs perm, \
  stat/IGAP/$(qtl_data)/$(chr).$(stat)_bed.gz \
  stat/IGAP/ld/$(qtl_data)/$(chr).$(stat)_ld.gz )))

jobs/step2-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/obs-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N med.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

jobs/step2-perm-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step2/perm-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N perm.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

# % = $(chr)
jobs/step2/obs-%.jobs: batches/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@for qtl_data in $(QTL_DATA) ; do zcat $< | awk -F'/' -vOFS=" " -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '{ batch=$$2; print "./make.qtl.assigned.R" OFS ("qtl/" chr "/" batch ".qtl-" qtl ".gz") OFS ("qtl/" chr "/" batch ".snps.gz") OFS ("qtl/" chr "/" batch ".genes.gz") OFS ("data/" gwas "/chr" chr ".ld_bed.gz") OFS ("$(TEMP)/matched/obs/" gwas "/" qtl "/" chr "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

# % = $(chr)
jobs/step2/perm-%.jobs: batches/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@for qtl_data in $(QTL_DATA) ; do zcat $< | awk -F'/' -vOFS=" " -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '{ batch=$$2; print "./make.qtl.assigned.R" OFS ("qtl/" chr "/" batch ".perm-" qtl ".gz") OFS ("qtl/" chr "/" batch ".snps.gz") OFS ("qtl/" chr "/" batch ".genes.gz") OFS ("data/" gwas "/chr" chr ".ld_bed.gz") OFS ("$(TEMP)/matched/perm/" gwas "/" qtl "/" chr "/" batch ".data.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

batches/%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 qtl/$*/b*.genes.gz | awk -F'/' '{ gsub(/.genes.gz/, "", $$NF); print $$(NF-1) "/" $$NF }' | gzip > $@

# % = $(qtl_data)/$(chr)
stat/IGAP/%.obs_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/obs/IGAP/$*/b*.data.gz | xargs zcat) | awk -F'\t' '{ printf $$1 FS int($$2) FS int($$3); for(j=4; j<=NF; ++j) printf FS $$j; printf "\n" }' | sort -k1,1 -k2,2n | gzip > $@

# % = $(qtl_data)/$(chr)
stat/IGAP/ld/%.obs_ld.gz : stat/IGAP/%.obs_bed.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk -F'\t' '{ ld[$$12] ++ } END { for(l in ld) print l FS ld[l] }' | sed 's/:/\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > $@

# % = $(qtl_data)/$(chr)
stat/IGAP/%.perm_bed.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	(ls -1 $(TEMP)/matched/perm/IGAP/$*/b*.data.gz | xargs zcat) | awk -F'\t' '{ printf $$1 FS int($$2) FS int($$3); for(j=4; j<=NF; ++j) printf FS $$j; printf "\n" }' | sort -k1,1 -k2,2n | gzip > $@

# % = $(qtl_data)/$(chr)
stat/IGAP/ld/%.perm_ld.gz : stat/IGAP/%.perm_bed.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk -F'\t' '{ ld[$$12] ++ } END { for(l in ld) print l FS ld[l] }' | sed 's/:/\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > $@


################################################################
# distribute LD jobs
step3: jobs/step3-jobs.txt.gz

step3-resubmit:
	zcat jobs/step3-jobs.txt.gz | awk '(system("[ ! -f " $$NF ".zqtl.txt.gz ]") == 0) && (system("[ ! -f " $$NF " ]") == 0)' | gzip > jobs/step3-jobs-resubmit.txt.gz
	@[ $$(zcat jobs/step3-jobs-resubmit.txt.gz | wc -l) -lt 1 ] || qsub -t 1-$$(zcat jobs/step3-jobs-resubmit.txt.gz | wc -l) -N ZQTL -binding "linear:1" -q long -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh jobs/step3-jobs-resubmit.txt.gz

jobs/step3-jobs.txt.gz: $(foreach chr, $(CHR), jobs/step3/obs-$(chr).jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N ZQTL -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	@rm $^

jobs/step3/obs-%.jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@

	@for qtl_data in $(QTL_DATA) ; do zcat stat/IGAP/ld/$${qtl_data}/$*.obs_ld.gz | awk -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.null.R" OFS ("stat/" gwas "/ld/" qtl "/" chr ".obs_ld.gz") OFS ("stat/" gwas "/" qtl "/" chr ".obs_bed.gz") FS ("1kg/chr" chr) FS NR FS 74000 FS 400 FS ("zqtl.null/IGAP/1kg/" qtl "/" chr "/" NR ".null.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

	@for qtl_data in $(QTL_DATA) ; do zcat stat/IGAP/ld/$${qtl_data}/$*.obs_ld.gz | awk -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.null.R" OFS ("stat/" gwas "/ld/" qtl "/" chr ".obs_ld.gz") OFS ("stat/" gwas "/" qtl "/" chr ".obs_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74000 FS 300 FS ("zqtl.null/IGAP/rosmap/" qtl "/" chr "/" NR ".null.gz") }' | awk 'system("[ ! -f " $$NF " ]") == 0' >> $@; done

	@for qtl_data in $(QTL_DATA) ; do zcat stat/IGAP/ld/$${qtl_data}/$*.obs_ld.gz | awk -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" gwas "/ld/" qtl "/" chr ".obs_ld.gz") OFS ("stat/" gwas "/" qtl "/" chr ".obs_bed.gz") FS ("1kg/chr" chr) FS NR FS 74000 FS 400 FS ("zqtl/IGAP/1kg/" qtl "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".zqtl.gz ]") == 0' >> $@; done

	@for qtl_data in $(QTL_DATA) ; do zcat stat/IGAP/ld/$${qtl_data}/$*.obs_ld.gz | awk -vchr=$* -vgwas=IGAP -vqtl=$${qtl_data} '$$(NF) >= 100 { print "./make.mediation.R" OFS ("stat/" gwas "/ld/" qtl "/" chr ".obs_ld.gz") OFS ("stat/" gwas "/" qtl "/" chr ".obs_bed.gz") FS ("geno/rosmap1709-chr" chr) FS NR FS 74000 FS 300 FS ("zqtl/IGAP/rosmap/" qtl "/" chr "/" NR) }' | awk 'system("[ ! -f " $$NF ".zqtl.gz ]") == 0' >> $@; done


################################################################
## Utilities
PLINKZIP := https://www.cog-genomics.org/static/bin/plink170906/plink_linux_x86_64.zip

bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/

