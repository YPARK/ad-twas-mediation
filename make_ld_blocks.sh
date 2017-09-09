#!/bin/bash -l

if [ $# -lt 4 ]; then
    exit 1
fi

chr=$1      # e.g., 1
ld=$2       # e.g., LD/igap-ad.bed.gz
gwas=$3     # e.g., ./IGAP/chr${chr}.txt.gz
outfile=$4

source /broad/software/scripts/useuse > /dev/null
# reuse -q R-3.3
# reuse -q GCC-5.2
# reuse -q .icc-2015
reuse -q BEDTools

# We computed pairwise correlations between pairs of variants in the
# Thousand Genomes European samples within 1 megabase and with r2 >
# 0.1.  We pruned to a desired threshold by iteratively picking the
# top-scoring variant (breaking ties arbitrarily) and removing the
# tagged variants until no variants remained.  (From Sarkar et al. 2016)

# We then merged the LD blocks if adjacent blocks overlap in genomic
# locations.
[ -f $outfile ] || \
    zcat $ld | \
	awk -F'\t' -vchr="chr"${chr} '$1 == chr {
  if(!($4 in lb) || (lb[$4] > $3)) lb[$4] = $3
  if(!($4 in ub) || (ub[$4] < $3)) ub[$4] = $3
  n[$4]++
  loc[$4] = (($4 in loc)? loc[$4] "|" : "") $3
}
END {
  for(k in lb) print chr FS lb[k] FS ub[k] FS n[k] FS loc[k]
}' | awk '$2 < $3' | sort -k1,1 -k2,2n | \
	bedtools merge -c 4,5 -o sum,distinct -delim "|" | \
	awk -F'\t' '{
  split($5,snps,"|")
  key = $1 ":" $2 ":" $3 ":" $4
  for(j=1; j<=$4; ++j)
    print $1 FS (snps[j] - 1) FS snps[j] FS key
}' | sort -k1,1 -k2,2n | \
	bedtools intersect -b stdin -a $gwas -wa -wb | \
	cut -f1-9,13 | sort -k1,1 -k2,2n | gzip > $outfile
