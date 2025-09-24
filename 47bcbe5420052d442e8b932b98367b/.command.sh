#!/bin/bash -ue
set -euo pipefail
samtools index aln.ont.bam
samtools faidx BIDI_mpox_case2_hybrid_final.fa
conda activate clair3
bash run_clair3.sh --bam_fn=aln.ont.bam --ref_fn=BIDI_mpox_case2_hybrid_final.fa --model_path="/opt/models/r1041_e82_400bps_sup_v430" --threads=32 --platform=ont --output=clair3_out
bcftools view -Oz -o ont.vcf.gz clair3_out/merge_output.vcf.gz
bcftools index ont.vcf.gz
bcftools norm -f BIDI_mpox_case2_hybrid_final.fa -m -both ont.vcf.gz \
  | bcftools filter -i 'QUAL>=20 && TYPE="INDEL"' \
  | bcftools sort -Oz -o ont.norm.indels.vcf.gz

# Check if VCF has variants, if not create empty indexed VCF
if [ ! -s ont.norm.indels.vcf.gz ] || [ $(bcftools view -H ont.norm.indels.vcf.gz | wc -l) -eq 0 ]; then
  # Create empty VCF with proper header
  bcftools view -h ont.vcf.gz | bcftools sort -Oz -o ont.norm.indels.vcf.gz
fi

# Try to create index, if it fails create empty placeholder
bcftools index -t ont.norm.indels.vcf.gz || {
  echo "Warning: bcftools index failed, creating empty .tbi file"
  touch ont.norm.indels.vcf.gz.tbi
}
