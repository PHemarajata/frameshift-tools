#!/bin/bash -ue
set -euo pipefail
freebayes -f BIDI_mpox_case2_hybrid_final.fa --min-alternate-fraction 0.2 --pooled-continuous aln.ilmn.bam > ilmn.vcf
bcftools norm -f BIDI_mpox_case2_hybrid_final.fa -m -both ilmn.vcf \
  | bcftools filter -i 'QUAL>=20 && TYPE="INDEL"' \
  | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz

# Check if VCF has variants, if not create empty indexed VCF
if [ ! -s ilmn.norm.indels.vcf.gz ] || [ $(bcftools view -H ilmn.norm.indels.vcf.gz | wc -l) -eq 0 ]; then
  # Create empty VCF with proper header
  bcftools view -h ilmn.vcf | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz
fi

# Try to create index, if it fails create empty placeholder
bcftools index ilmn.norm.indels.vcf.gz || {
  echo "Warning: bcftools index failed, creating empty .tbi file"
  touch ilmn.norm.indels.vcf.gz.tbi
}
