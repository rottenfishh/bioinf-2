"""Сборка итоговой таблицы: VEP (vcf.gz) + ClinVar (tsv) + gnomAD (из VEP или gnomad_manual.tsv)."""
import gzip, csv, glob, os, sys
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# порядок тяжести последствий — по Ensembl VEP (от самого тяжёлого)
SEV = """transcript_ablation splice_acceptor_variant splice_donor_variant stop_gained frameshift_variant
stop_lost start_lost transcript_amplification feature_elongation feature_truncation inframe_insertion
inframe_deletion missense_variant protein_altering_variant splice_donor_5th_base_variant
splice_region_variant splice_donor_region_variant splice_polypyrimidine_tract_variant
incomplete_terminal_codon_variant start_retained_variant stop_retained_variant synonymous_variant
coding_sequence_variant mature_miRNA_variant 5_prime_UTR_variant 3_prime_UTR_variant
non_coding_transcript_exon_variant intron_variant NMD_transcript_variant non_coding_transcript_variant
coding_transcript_variant upstream_gene_variant downstream_gene_variant TFBS_ablation
TFBS_amplification TF_binding_site_variant regulatory_region_ablation regulatory_region_amplification
regulatory_region_variant intergenic_variant sequence_variant""".split()
rank = lambda cons: min(SEV.index(c) for c in cons.split("&"))

# ---------- VEP ----------
vep_file = sorted(glob.glob("student_37_*VEP*.vcf.gz"))[-1]  # берём последний запуск (с gnomAD)
fmt, variants, csq = None, [], {}
for line in gzip.open(vep_file, "rt"):
    if line.startswith("##INFO=<ID=CSQ"):
        fmt = line.split("Format: ")[1].split('"')[0].split("|")
    if line.startswith("#"):
        continue
    c, p, i, r, a, _, _, info = line.rstrip("\n").split("\t")[:8]
    key = (c, p, r, a); variants.append((c, p, i, r, a)); csq[key] = []
    s = [x for x in info.split(";") if x.startswith("CSQ=")][0][4:]
    for t in s.split(","):
        csq[key].append(dict(zip(fmt, t.split("|"))))

# длинная таблица: все транскрипты
with open("results/vep_all_transcripts.tsv", "w", newline="") as f:
    cols = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature","BIOTYPE","EXON","INTRON",
            "Protein_position","Amino_acids","Codons","STRAND","CANONICAL","MANE","GENCODE_PRIMARY"]
    extra = [x for x in ("gnomAD_exomes_AF","gnomAD_genomes_AF","gnomAD_genomes_AF_grpmax") if x in fmt]
    w = csv.writer(f, delimiter="\t")
    w.writerow(["CHROM","POS","ID","REF","ALT"] + cols + [x for x in extra if x not in cols])
    for c, p, i, r, a in variants:
        for d in csq[(c, p, r, a)]:
            w.writerow([c, p, i, r, a] + [d.get(x, "") for x in cols] + [d.get(x, "") for x in extra if x not in cols])

# ---------- ClinVar ----------
clin = {}
for f in glob.glob("clinvar_download_*.tsv"):
    for row in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"):
        _, pos0, ref, alt = row["Canonical SPDI"].split(":")
        key = (row["GRCh38Chromosome"], str(int(pos0) + 1), ref, alt)
        clin[key] = row

# ---------- gnomAD (ручная таблица, если есть) ----------

def af_from_vep(ts, field):
    vals = {d.get(field, "") for d in ts} - {""}
    return ";".join(sorted(vals))

def prot(d):
    if d["Amino_acids"] and d["Protein_position"]:
        aa = d["Amino_acids"].split("/")
        return f"p.{aa[0]}{d['Protein_position']}{aa[-1] if len(aa) > 1 else '='}"
    return ""

out = []
for n, (c, p, i, r, a) in enumerate(variants, 1):
    ts = csq[(c, p, r, a)]
    mane = [d for d in ts if d["MANE"]]
    pool, rule = (mane, "MANE Select") if mane else ([d for d in ts if d["CANONICAL"] == "YES"] or ts, "canonical")
    best = min(pool, key=lambda d: rank(d["Consequence"]))
    worst = min(ts, key=lambda d: rank(d["Consequence"]))
    cv = clin.get((c, p, r, a))
    vid = f"{c}-{p}-{r}-{a}"
    out.append({
        "N": n, "variant": f"{c}:{p} {r}>{a}", "rsID": i,
        "gene": best["SYMBOL"],
        "consequence_VEP": best["Consequence"], "impact": best["IMPACT"],
        "transcript": best["Feature"], "transcript_rule": rule, "protein_change": prot(best),
        "most_severe_any_transcript": worst["Consequence"],
        "most_severe_gene_transcript": f"{worst['SYMBOL']} {worst['Feature']} ({worst['BIOTYPE']})",
        "n_transcripts": len(ts),
        "gnomAD_exomes_AF": af_from_vep(ts, "gnomAD_exomes_AF") or "not found",
        "gnomAD_genomes_AF": af_from_vep(ts, "gnomAD_genomes_AF") or "not found",
        "gnomAD_genomes_AF_grpmax": af_from_vep(ts, "gnomAD_genomes_AF_grpmax"),
        "ClinVar": "yes" if cv else "no",
        "ClinVar_classification": cv["Germline classification"] if cv else "",
        "ClinVar_review_status": cv["Germline review status"] if cv else "",
        "ClinVar_conditions": cv["Condition(s)"] if cv else "",
        "ClinVar_accession": cv["Accession"] if cv else "",
    })

with open("results/annotated_variants.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
print(f"VEP: {vep_file}; ClinVar записей сопоставлено: {sum(o['ClinVar']=='yes' for o in out)} из {len(clin)}")
