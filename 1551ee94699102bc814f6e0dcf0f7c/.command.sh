#!/bin/bash -ue
set -euo pipefail
      export ALIAS_MAP=""

      # Intersections
      bedtools intersect -wa -wb -a ilmn.indels.bed -b cds.bed > ilmn.cds.tsv || true
      bedtools intersect -wa -wb -a ont.indels.bed  -b cds.bed > ont.cds.tsv  || true

      python3 - <<'PY'
import csv, json, os, re

# Load optional alias map (TSV: ncbi_name  nextclade_name)
alias_map = {}
alias_path = os.environ.get("ALIAS_MAP")
if alias_path and os.path.exists(alias_path):
    with open(alias_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            a=line.rstrip().split('\t')
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
            g = re.sub(r'^"|"$','',c).strip()
            return alias_map.get(g, g)
    return ""

def load_intersections(path):
    out=[]
    if not os.path.exists(path) or os.path.getsize(path)==0: return out
    with open(path) as f:
        for line in f:
            a=line.rstrip().split('\t')
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
