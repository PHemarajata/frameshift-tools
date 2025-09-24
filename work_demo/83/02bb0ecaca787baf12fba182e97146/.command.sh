#!/bin/bash -ue
set -euo pipefail
minimap2 -t 32 -ax map-ont BIDI_mpox_case2_hybrid_final.fa SImpxv_case2_nanopore.fastq.gz | samtools sort -@8 -o aln.ont.bam
samtools index aln.ont.bam
