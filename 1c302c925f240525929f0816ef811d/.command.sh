#!/bin/bash -ue
set -euo pipefail
: > ilmn.indels.bed; : > ont.indels.bed
if [ -f "ilmn.norm.indels.vcf.gz" ] && [ "ilmn.norm.indels.vcf.gz" != "NO_FILE" ]; then
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%QUAL\t[%AF]\t%TYPE\n' ilmn.norm.indels.vcf.gz \
    | awk 'BEGIN{OFS="\t"}{lr=length($3);la=length($4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print $1,$2-1,$2,"ILMN",$3,$4,d,fs,$5,$6,$7,$8}' > ilmn.indels.bed
fi
if [ -f "ont.norm.indels.vcf.gz" ] && [ "ont.norm.indels.vcf.gz" != "NO_FILE" ]; then
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%QUAL\t[%AF]\t%TYPE\n' ont.norm.indels.vcf.gz \
    | awk 'BEGIN{OFS="\t"}{lr=length($3);la=length($4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print $1,$2-1,$2,"ONT",$3,$4,d,fs,$5,$6,$7,$8}' > ont.indels.bed
fi
