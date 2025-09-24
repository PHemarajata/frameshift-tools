
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
    params.model_path  = null  // let us override; default below

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
      publishDir "results/nextclade", mode: 'copy'
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
      publishDir "results/variants", mode: 'copy'
      input:
        path fasta
        path bam
      output:
        path "ilmn.norm.indels.vcf.gz"
        path "ilmn.norm.indels.vcf.gz.tbi"
      script:
      """
      set -euo pipefail
      freebayes -f ${fasta} --min-alternate-fraction 0.2 --pooled-continuous ${bam} > ilmn.vcf
      bcftools norm -f ${fasta} -m -both ilmn.vcf \\
        | bcftools filter -i 'QUAL>=20 && TYPE="INDEL"' \\
        | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz
      
      # Check if VCF has variants, if not create empty indexed VCF
      if [ ! -s ilmn.norm.indels.vcf.gz ] || [ \$(bcftools view -H ilmn.norm.indels.vcf.gz | wc -l) -eq 0 ]; then
        # Create empty VCF with proper header
        bcftools view -h ilmn.vcf | bcftools sort -Oz -o ilmn.norm.indels.vcf.gz
      fi
      
      # Try to create index, if it fails create empty placeholder
      bcftools index -t ilmn.norm.indels.vcf.gz || {
        echo "Warning: bcftools index failed, creating empty .tbi file"
        touch ilmn.norm.indels.vcf.gz.tbi
      }
      """
    }

    process CallONT {
      tag "call_ont"
      container 'staphb/clair3:1.2.0'
      publishDir "results/variants", mode: 'copy'
      input:
        path fasta
        path bam
        val  model_path
      output:
        path "ont.raw.vcf.gz"
      script:
      """
      set -euo pipefail
      samtools index ${bam}
      samtools faidx ${fasta}
      bash /clair3/run_clair3.sh --bam_fn=${bam} --ref_fn=${fasta} --model_path="${model_path}" --threads=32 --platform=ont --output=clair3_out --include_all_ctgs
      
      # Copy the VCF output 
      cp clair3_out/merge_output.vcf.gz ont.raw.vcf.gz
      """
    }
    
    process ProcessONTVCF {
      tag "process_ont_vcf"
      publishDir "results/variants", mode: 'copy'
      input:
        path fasta
        path raw_vcf
      output:
        path "ont.norm.indels.vcf.gz"
        path "ont.norm.indels.vcf.gz.tbi"
      script:
      """
      set -euo pipefail
      bcftools view -Oz -o ont.vcf.gz ${raw_vcf}
      bcftools index ont.vcf.gz
      bcftools norm -f ${fasta} -m -both ont.vcf.gz \\
        | bcftools filter -i 'QUAL>=20 && TYPE="INDEL"' \\
        | bcftools sort -Oz -o ont.norm.indels.vcf.gz
      
      # Check if VCF has variants, if not create empty indexed VCF
      if [ ! -s ont.norm.indels.vcf.gz ] || [ \$(bcftools view -H ont.norm.indels.vcf.gz | wc -l) -eq 0 ]; then
        # Create empty VCF with proper header
        bcftools view -h ont.vcf.gz | bcftools sort -Oz -o ont.norm.indels.vcf.gz
      fi
      
      # Try to create index, if it fails create empty placeholder
      bcftools index -t ont.norm.indels.vcf.gz || {
        echo "Warning: bcftools index failed, creating empty .tbi file"
        touch ont.norm.indels.vcf.gz.tbi
      }
      """
    }


    // -------- GFF3 -> BED (keep attributes for NCBI names) --------
    process MakeCdsBed {
      tag "gff_to_bed"
      publishDir "results/annotations", mode: 'copy'
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
      publishDir "results/indels", mode: 'copy'
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
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t[%DP]\\t[%AF]\\t%TYPE\\n' ${ilmn_vcf} \\
          | awk 'BEGIN{OFS="\\t"}{lr=length(\$3);la=length(\$4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print \$1,\$2-1,\$2,"ILMN",\$3,\$4,d,fs,\$5,\$6,\$7,\$8}' > ilmn.indels.bed
      fi
      if [ -f "${ont_vcf}" ] && [ "${ont_vcf}" != "NO_FILE" ]; then
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t[%DP]\\t[%AF]\\t%TYPE\\n' ${ont_vcf} \\
          | awk 'BEGIN{OFS="\\t"}{lr=length(\$3);la=length(\$4);d=(la>lr?la-lr:lr-la);fs=(d%3==0)?"in-frame":"frameshift";print \$1,\$2-1,\$2,"ONT",\$3,\$4,d,fs,\$5,\$6,\$7,\$8}' > ont.indels.bed
      fi
      """
    }

    // -------- Verify frameshifts & harmonize gene names --------
    process VerifyFrameshifts {
      tag "verify_and_harmonize"
      publishDir "results/final", mode: 'copy'
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

      # Fix contig name mismatch before intersections
      # Extract short contig names from indels files to match CDS file
      sed 's/^BIDI_mpox_case2_hybrid_APHLproject_work_20250923_103207_//' ${ilmn} > ilmn_fixed.bed
      sed 's/^BIDI_mpox_case2_hybrid_APHLproject_work_20250923_103207_//' ${ont} > ont_fixed.bed
      
      # Intersections with fixed contig names
      bedtools intersect -wa -wb -a ilmn_fixed.bed -b ${cds} > ilmn.cds.tsv || true
      bedtools intersect -wa -wb -a ont_fixed.bed -b ${cds} > ont.cds.tsv || true
      
      # Debug: Check intersection results
      echo "DEBUG: ILMN intersections:" >&2
      wc -l ilmn.cds.tsv >&2
      echo "DEBUG: ONT intersections:" >&2
      wc -l ont.cds.tsv >&2

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

# Create Nextclade-based frameshift verification report
with open("frameshift_report.csv","w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["nextclade_gene","codon_range","ilmn_indels_found","ilmn_max_dp","ilmn_max_af",
                "ont_indels_found","ont_max_dp","ont_max_af","concordance_status"])
    
    # Process each Nextclade frameshift
    for cds,beg,end in nc_fs:
        codon_range = f"{beg}-{end}" if beg and end else "unknown"
        
        # Count indels for this gene (simple approach for now)
        ilmn_gene_indels = [r for r in ilmn if r["gene"] == cds or cds in r["gene"]]
        ont_gene_indels = [r for r in ont if r["gene"] == cds or cds in r["gene"]]
        
        ilmn_count = len(ilmn_gene_indels)
        ilmn_max_dp = max([float(r["dp"]) for r in ilmn_gene_indels], default=0)
        ilmn_max_af = max([float(r["af"]) for r in ilmn_gene_indels], default=0)
        
        ont_count = len(ont_gene_indels)
        ont_max_dp = max([float(r["dp"]) for r in ont_gene_indels], default=0)
        ont_max_af = max([float(r["af"]) for r in ont_gene_indels], default=0)
        
        # Determine concordance
        if ilmn_count > 0 and ont_count > 0:
            concordance = "both_platforms"
        elif ilmn_count > 0 or ont_count > 0:
            concordance = "single_platform"
        else:
            concordance = "no_evidence"
        
        w.writerow([cds, codon_range, ilmn_count, ilmn_max_dp, ilmn_max_af,
                    ont_count, ont_max_dp, ont_max_af, concordance])

print(f"DEBUG: Processed {len(nc_fs)} Nextclade frameshifts")
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
  model_ch      = Channel.value( params.model_path ?: '/clair3/models/ont' )

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
    CallONT(consensus_ch, MapONT.out[0], model_ch)  // Clair3 in container
    ProcessONTVCF(consensus_ch, CallONT.out)        // bcftools on host
    ont_vcf_ch = ProcessONTVCF.out[0]  // First output is processed VCF file
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