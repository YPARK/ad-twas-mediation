
CHR := $(shell seq 1 22)

step1: $(foreach chr, $(CHR), data/IGAP/chr$(chr).ld_bed.gz)

# % = $(chr)
data/IGAP/chr%.ld_bed.gz: LD/igap-ad.bed.gz IGAP/chr%.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make_ld_blocks.sh $* $^ $@

