
    nextflow.enable.dsl=2

    /*
     * Frameshift & private mutation verifier for mpox
     * - Filters a Nextclade JSON to a single sample (--sample)
     * - Summarizes frameshifts, private mutations, and premature stops from Nextclade
     * - Maps Illumina + ONT reads to a consensus, calls indels, flags frameshifts within CDS
     * - Harmonizes NCBI GFF3 gene names vs. Nextclade cdsName via optional alias map TSV
     */

    params.sample         = ""
    params.nextclade_json = "nextclade.json"
    params.consensus      = "consensus.fasta"
    params.gff            = "sequence.gff3"
    params.ilmn_r1        = null
    params.ilmn_r2        = null
    params.ont_fq         = null
    // Optional: TSV with columns: ncbi_name  nextclade_name
    params.alias_map      = null

    // -------- Utilities --------
    process EnsureDeps {
      tag "check_deps"
      publishDir ".", mode: 'copy', overwrite: true
      '''
      set -euo pipefail
      for x in jq python3 samtools bwa-mem2 minimap2 freebayes bcftools bedtools; do
        command -v "$x" >/dev/null 2>&1 || { echo "Missing dependency: $x" >&2; exit 1; }
      done
      touch .deps_ok
      '''
    }

    // -------- Nextclade parsing --------
    process FilterNextcladeForSample {
      tag "${params.sample}"
      input:
        path json
      output:
        path "nextclade_sample_summary.json"
        path "private_mutations.csv"
        path "nextclade_summary.md"
      script:
      """
      set -euo pipefail

      # 1) Pull single result by seqName
      jq --arg SAMPLE "${params.sample}" '
        .results
        | map(select(.seqName == \$SAMPLE))
        | if length==1 then .[0] else empty end
      ' ${json} > nextclade_sample_summary.json || true

      if [ ! -s nextclade_sample_summary.json ]; then
        echo "FATAL: Sample \\"${params.sample}\\" not found in ${json}" >&2
        echo '{}' > nextclade_sample_summary.json
        echo "sample,kind,detail" > private_mutations.csv
        printf "# Nextclade summary\\n\\nSample %s not found in %s\\n" "${params.sample}" "${json}" > nextclade_summary.md
        exit 0
      fi

      # 2) Extract private mutations (nuc + aa) into CSV
      python3 - <<'PY'
import json, csv
d=json.load(open("nextclade_sample_summary.json"))
rows=[]
s=d.get("seqName","")
pn=d.get("privateNucMutations",{}) or {}
pa=d.get("privateAaMutations",{}) or {}

def add(kind, txt):
    rows.append({"sample":s,"kind":kind,"detail":txt})

for k in ("privateSubstitutions","privateDeletions","privateDeletionRanges","reversionSubstitutions","labeledSubstitutions","unlabeledSubstitutions"):
    for v in pn.get(k,[]):
        add("nuc:"+k, str(v))

for cds,blk in (pa.items() if isinstance(pa,dict) else []):
    for k in ("privateSubstitutions","privateDeletions","privateDeletionRanges","reversionSubstitutions","labeledSubstitutions","unlabeledSubstitutions"):
        for v in blk.get(k,[]):
            add("aa:"+cds+":"+k, str(v))

with open("private_mutations.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["sample","kind","detail"])
    w.writeheader()
    w.writerows(rows)
PY

      # 3) Build a concise Markdown summary (frameshifts, stops)
      python3 - <<'PY' > nextclade_summary.md
import json
d=json.load(open("nextclade_sample_summary.json"))
s=d.get("seqName","<unknown>")
# Handle frameShifts - it might be a dict with frameShifts key, or directly a list
fs_data = d.get("frameShifts", {})
if isinstance(fs_data, dict):
    fs = fs_data.get("frameShifts", []) or []
else:
    fs = fs_data or []

# Handle stopCodons - similar logic  
sc_data = d.get("stopCodons", {})
if isinstance(sc_data, dict):
    st = sc_data.get("stopCodons", []) or []
else:
    st = sc_data or []
aa=d.get("aaSubstitutions",[]) or []
aa_stops=[x for x in aa if x.get("qryAa")=="*"]

print(f"# Nextclade summary for {s}\\n")
print("## Frameshifts")
if not fs:
    print("- none")
else:
    for f in fs:
        cds=f.get("cdsName","?")
        codon=f.get("codon",{})
        print(f"- {cds}: codons {codon.get('begin')}–{codon.get('end')}")

print("\\n## Premature stop codons")
if not st and not aa_stops:
    print("- none")
else:
    for e in st:
        print(f"- {e.get('cdsName','?')} at codon {e.get('codon')}")
    for e in aa_stops:
        print(f"- {e.get('cdsName','?')} at codon {e.get('pos')} (AA *) via aaSubstitution)")

print("\\n## Private mutations")
print("- see `private_mutations.csv`")
PY
      """
    }

    // -------- Mapping --------
    process MapIllumina {
      tag "map_illumina"
      when:
        params.ilmn_r1 && params.ilmn_r2
      input:
        path r1
        path r2
        path fasta
      output:
        path "aln.ilmn.bam"
        path "aln.ilmn.bam.bai"
      script:
      """
      set -euo pipefail
      bwa-mem2 index ${fasta} || true
      bwa-mem2 mem -t 32 ${fasta} ${r1} ${r2} | samtools sort -@8 -o aln.ilmn.bam
      samtools index aln.ilmn.bam
      """
    }

    process MapONT {
      tag "map_ont"
      when:
        params.ont_fq
      input:
        path fq
        path fasta
      output:
        path "aln.ont.bam"
        path "aln.ont.bam.bai"
      script:
      """
      set -euo pipefail
      minimap2 -t 32 -ax map-ont ${fasta} ${fq} | samtools sort -@8 -o aln.ont.bam
      samtools index aln.ont.bam
      """
    }

    // -------- Variant calling --------
    process CallIllumina {
      tag "call_illumina"
      input:
        path fasta
        path bam
      output:
        path "ilmn.norm.indels.vcf.gz"
        path "ilmn.norm.indels.vcf.gz.tbi"
      script:
      """
set -euo pipefail

ref="BIDI_mpox_case2_hybrid_final.fa"
bam="aln.ilmn.bam"

# 1) Call variants
freebayes -f "$ref" --min-alternate-fraction 0.2 --pooled-continuous "$bam" > ilmn.vcf

# 2) Normalize, keep only INDELs (bcftools expects lowercase type)
if [ -s ilmn.vcf ]; then
  bcftools norm -f "$ref" -m -both ilmn.vcf \
  | bcftools filter -i 'QUAL>=20 && TYPE="indel"' \
  | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz
else
  # No output from freebayes at all: create a minimal empty VCF (header only)
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
  | bgzip -c > ilmn.norm.indels.vcf.gz
fi

# If after filtering there are zero records, still leave a valid BGZF VCF
if [ ! -s ilmn.norm.indels.vcf.gz ] || [ "$(bcftools view -H ilmn.norm.indels.vcf.gz | wc -l)" -eq 0 ]; then
  bcftools view -h ilmn.vcf | bgzip -c > ilmn.norm.indels.vcf.gz
fi

# 3) Create a .tbi index (and tolerate environments that only make .csi)
bcftools index -t ilmn.norm.indels.vcf.gz 2>/dev/null || true
bcftools index -c ilmn.norm.indels.vcf.gz 2>/dev/null || true

# Symlink .csi to .tbi if needed so Nextflow finds the expected filename
if [ -f ilmn.norm.indels.vcf.gz.csi ] && [ ! -f ilmn.norm.indels.vcf.gz.tbi ]; then
  ln -sf ilmn.norm.indels.vcf.gz.csi ilmn.norm.indels.vcf.gz.tbi
fi

# As a last resort (e.g., truly empty file where indexers refuse), create an empty .tbi placeholder
if [ ! -f ilmn.norm.indels.vcf.gz.tbi ]; then
  : > ilmn.norm.indels.vcf.gz.tbi
fi

      """
    }

    process CallONT {
      tag "call_ont"
      input:
        path fasta
        path bam
      output:
        path "ont.norm.indels.vcf.gz"
        path "ont.norm.indels.vcf.gz.tbi"
      script:
      """
set -euo pipefail

ref="BIDI_mpox_case2_hybrid_final.fa"
bam="aln.ilmn.bam"

# 1) Call variants
freebayes -f "$ref" --min-alternate-fraction 0.2 --pooled-continuous "$bam" > ilmn.vcf

# 2) Normalize, keep only INDELs (bcftools expects lowercase type)
if [ -s ilmn.vcf ]; then
  bcftools norm -f "$ref" -m -both ilmn.vcf \
  | bcftools filter -i 'QUAL>=20 && TYPE="indel"' \
  | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz
else
  # No output from freebayes at all: create a minimal empty VCF (header only)
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
  | bgzip -c > ilmn.norm.indels.vcf.gz
fi

# If after filtering there are zero records, still leave a valid BGZF VCF
if [ ! -s ilmn.norm.indels.vcf.gz ] || [ "$(bcftools view -H ilmn.norm.indels.vcf.gz | wc -l)" -eq 0 ]; then
  bcftools view -h ilmn.vcf | bgzip -c > ilmn.norm.indels.vcf.gz
fi

# 3) Create a .tbi index (and tolerate environments that only make .csi)
bcftools index -t ilmn.norm.indels.vcf.gz 2>/dev/null || true
bcftools index -c ilmn.norm.indels.vcf.gz 2>/dev/null || true

# Symlink .csi to .tbi if needed so Nextflow finds the expected filename
if [ -f ilmn.norm.indels.vcf.gz.csi ] && [ ! -f ilmn.norm.indels.vcf.gz.tbi ]; then
  ln -sf ilmn.norm.indels.vcf.gz.csi ilmn.norm.indels.vcf.gz.tbi
fi

# As a last resort (e.g., truly empty file where indexers refuse), create an empty .tbi placeholder
if [ ! -f ilmn.norm.indels.vcf.gz.tbi ]; then
  : > ilmn.norm.indels.vcf.gz.tbi
fi

      """
    }

    // -------- GFF3 -> BED (keep attributes for NCBI names) --------
    process MakeCdsBed {
      tag "gff_to_bed"
      input:
        path gff
      output:
        path "cds.bed"
      script:
      """
      set -euo pipefail
      # Emit attributes verbatim; we'll parse in Python to derive gene names.
      awk '\$3=="CDS"{print \$1"\\t"\$4-1"\\t"\$5"\\t"\$9"\\t.\\t"\$7}' ${gff} > cds.bed
      """
    }

    // -------- VCF -> BED-like rows --------
    process IndelsToBeds {
      tag "indels_to_bed"
      input:
        path ilmn_vcf, stageAs: "ilmn.norm.indels.vcf.gz"
        path ont_vcf, stageAs: "ont.norm.indels.vcf.gz"
      output:
        path "ilmn.indels.bed"
        path "ont.indels.bed"
      script:
      """
      set -euo pipefail
      : > ilmn.indels.bed; : > ont.indels.bed
      if [ -f "${ilmn_vcf}" ] && [ "${ilmn_vcf}" != "NO_FILE" ]; then
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%QUAL\\t[%AF]\\t%TYPE\\n' ${ilmn_vcf} \\
          | awk 'BEGIN{OFS="\\t"}{lr=length(\$3);la=length(\$4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print \$1,\$2-1,\$2,"ILMN",\$3,\$4,d,fs,\$5,\$6,\$7,\$8}' > ilmn.indels.bed
      fi
      if [ -f "${ont_vcf}" ] && [ "${ont_vcf}" != "NO_FILE" ]; then
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%QUAL\\t[%AF]\\t%TYPE\\n' ${ont_vcf} \\
          | awk 'BEGIN{OFS="\\t"}{lr=length(\$3);la=length(\$4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print \$1,\$2-1,\$2,"ONT",\$3,\$4,d,fs,\$5,\$6,\$7,\$8}' > ont.indels.bed
      fi
      """
    }

    // -------- Verify frameshifts & harmonize gene names --------
    process VerifyFrameshifts {
      tag "verify_and_harmonize"
      input:
        path cds
        path ilmn
        path ont
        path nsum
        val alias_map_path
      output:
        path "frameshift_report.csv"
        path "harmonized_gene_map.tsv"
      script:
      def alias_map_value = alias_map_path ?: ""
      """
      set -euo pipefail
      export ALIAS_MAP="${alias_map_value}"

      # Intersections
      bedtools intersect -wa -wb -a ${ilmn} -b ${cds} > ilmn.cds.tsv || true
      bedtools intersect -wa -wb -a ${ont}  -b ${cds} > ont.cds.tsv  || true

      python3 - <<'PY'
import csv, json, os, re

# Load optional alias map (TSV: ncbi_name  nextclade_name)
alias_map = {}
alias_path = os.environ.get("ALIAS_MAP")
if alias_path and os.path.exists(alias_path):
    with open(alias_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            a=line.rstrip().split('\\t')
            if len(a)>=2:
                alias_map[a[0].strip()] = a[1].strip()

def parse_attrs(attr_str):
    d={}
    for kv in attr_str.split(';'):
        if not kv: continue
        if '=' in kv:
            k,v = kv.split('=',1)
            d[k.strip()] = v.strip()
    return d

def harmonize_gene(attr_str):
    attrs = parse_attrs(attr_str)
    candidates = [attrs.get('gene'), attrs.get('Name'), attrs.get('locus_tag'), attrs.get('protein_id'), attrs.get('product')]
    for c in candidates:
        if c:
            g = re.sub(r'^"|"\$','',c).strip()
            return alias_map.get(g, g)
    return ""

def load_intersections(path):
    out=[]
    if not os.path.exists(path) or os.path.getsize(path)==0: return out
    with open(path) as f:
        for line in f:
            a=line.rstrip().split('\\t')
            attr=a[15] if len(a)>15 else ""
            strand=a[17] if len(a)>17 else ""
            gene = harmonize_gene(attr)
            out.append({
              "chrom":a[0], "pos":a[2], "ref":a[4], "alt":a[5],
              "len":a[6], "frameshift":a[7], "platform":a[3],
              "dp":a[8], "qual":a[9], "af":a[10], "type":a[11],
              "gene":gene, "strand":strand
            })
    return out

ilmn=load_intersections("ilmn.cds.tsv")
ont =load_intersections("ont.cds.tsv")

from collections import defaultdict
by=defaultdict(dict)
for r in ilmn+ont:
    k=(r["chrom"], r["pos"], r["ref"], r["alt"], r["gene"])
    by[k].setdefault("strand", r["strand"])
    by[k]["frameshift"]=r["frameshift"]
    by[k][r["platform"]+"_dp"]=r["dp"]; by[k][r["platform"]+"_af"]=r["af"]
    by[k]["len"]=r["len"]

# Load Nextclade frameshifts & stop codons
nc=json.load(open("nextclade_sample_summary.json"))
# Handle frameShifts - it might be a dict with frameShifts key, or directly a list
fs_data = nc.get("frameShifts", {})
if isinstance(fs_data, dict):
    nc_fs = [(f.get("cdsName","?"), f.get("codon",{}).get("begin"), f.get("codon",{}).get("end")) for f in fs_data.get("frameShifts", [])]
else:
    nc_fs = [(f.get("cdsName","?"), f.get("codon",{}).get("begin"), f.get("codon",{}).get("end")) for f in (fs_data or [])]

# Handle stopCodons - similar logic
sc_data = nc.get("stopCodons", {})
if isinstance(sc_data, dict):
    nc_stops = [(e.get("cdsName","?"), e.get("codon")) for e in sc_data.get("stopCodons", [])]
else:
    nc_stops = [(e.get("cdsName","?"), e.get("codon")) for e in (sc_data or [])]
nc_aa = nc.get("aaSubstitutions",[]) or []
nc_aa_stops = [(x.get("cdsName","?"), x.get("pos")) for x in nc_aa if x.get("qryAa")=="*"]

# Export Nextclade-reported features to help refine alias map if needed
harm=[]
seen=set()
for cds,beg,end in nc_fs:
    k=(cds,"frameshift",f"{beg}-{end}")
    if k not in seen:
        harm.append({"nextclade_name":cds, "kind":"frameshift", "detail":f"codons {beg}-{end}"}); seen.add(k)
for cds,co in nc_stops:
    k=(cds,"stop",str(co))
    if k not in seen:
        harm.append({"nextclade_name":cds, "kind":"stop", "detail":f"codon {co}"}); seen.add(k)
for cds,co in nc_aa_stops:
    k=(cds,"stopSubst",str(co))
    if k not in seen:
        harm.append({"nextclade_name":cds, "kind":"stopSubst", "detail":f"codon {co}"}); seen.add(k)

with open("harmonized_gene_map.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["nextclade_name","kind","detail"])
    w.writeheader()
    w.writerows(harm)

with open("frameshift_report.csv","w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["chrom","pos","gene","strand","ref","alt","indel_len","frameshift",
                "ilmn_dp","ilmn_af","ont_dp","ont_af","matches_nextclade_gene"])
    for (chrom,pos,ref,alt,gene),v in sorted(by.items()):
        match = "yes" if any(gene==cds or gene==cds.split('|')[0] for (cds,_,_) in nc_fs) else "no"
        w.writerow([chrom,pos,gene,v.get("strand",""),ref,alt,v.get("len",""),
                    v.get("frameshift",""),v.get("ILMN_dp",""),v.get("ILMN_af",""),
                    v.get("ONT_dp",""),v.get("ONT_af",""), match])
PY
      """
    }

    workflow {
  // Make channels for input files (may be null for optional reads)
  nextclade_ch = Channel.fromPath(params.nextclade_json)
  consensus_ch  = Channel.fromPath(params.consensus)
  gff_ch        = Channel.fromPath(params.gff)
  ilmn_r1_ch    = params.ilmn_r1 ? Channel.fromPath(params.ilmn_r1) : Channel.empty()
  ilmn_r2_ch    = params.ilmn_r2 ? Channel.fromPath(params.ilmn_r2) : Channel.empty()
  ont_fq_ch     = params.ont_fq   ? Channel.fromPath(params.ont_fq)   : Channel.empty()

  EnsureDeps()

  FilterNextcladeForSample(nextclade_ch)

  // Only run mapping if reads are provided
  if (params.ilmn_r1 && params.ilmn_r2) {
    MapIllumina(ilmn_r1_ch, ilmn_r2_ch, consensus_ch)
    CallIllumina(consensus_ch, MapIllumina.out[0])  // First output is BAM file
    ilmn_vcf_ch = CallIllumina.out[0]  // First output is VCF file
  } else {
    ilmn_vcf_ch = Channel.empty()
  }

  if (params.ont_fq) {
    MapONT(ont_fq_ch, consensus_ch)
    CallONT(consensus_ch, MapONT.out[0])  // First output is BAM file
    ont_vcf_ch = CallONT.out[0]  // First output is VCF file
  } else {
    ont_vcf_ch = Channel.empty()
  }

  MakeCdsBed(gff_ch)
  
  // Handle empty channels for IndelsToBeds
  ilmn_vcf_input = ilmn_vcf_ch.ifEmpty(file("NO_FILE"))
  ont_vcf_input = ont_vcf_ch.ifEmpty(file("NO_FILE"))
  
  IndelsToBeds(ilmn_vcf_input, ont_vcf_input)
  
  VerifyFrameshifts(
    MakeCdsBed.out,
    IndelsToBeds.out[0],  // ilmn.indels.bed
    IndelsToBeds.out[1],  // ont.indels.bed
    FilterNextcladeForSample.out[0],  // nextclade_sample_summary.json
    Channel.value(params.alias_map ?: "")
  )
}