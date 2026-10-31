import gzip, csv, sys, glob
from collections import Counter
f = glob.glob("student_37_*_VEP.vcf.gz")[0]
fmt=None; rows=[]
for line in gzip.open(f,"rt"):
    if line.startswith("##INFO=<ID=CSQ"):
        fmt=line.split("Format: ")[1].split('"')[0].split("|")
    if line.startswith("#"): continue
    c,p,i,r,a,_,_,info=line.rstrip("\n").split("\t")[:8]
    csq=[x for x in info.split(";") if x.startswith("CSQ=")][0][4:]
    for t in csq.split(","):
        d=dict(zip(fmt,t.split("|")))
        d.update(CHROM=c,POS=p,ID=i,REF=r,ALT=a); rows.append(d)
print(fmt); print(len(rows))
print(Counter(d["Feature_type"] for d in rows), Counter(d["SOURCE"] for d in rows))
print(Counter(c for d in rows for c in d["Consequence"].split("&")))
ids=Counter(d["ID"] for d in rows)
mane={d["ID"] for d in rows if d["MANE_SELECT"]}
print("no MANE:", [i for i in ids if i not in mane])
print("no MANE:", [i for i in ids if i not in {d["ID"] for d in rows if d["MANE"]}])
for d in rows:
    if "stop_gained" in d["Consequence"] or "splice" in d["Consequence"]:
        print(d["ID"],d["SYMBOL"],d["Feature"],d["BIOTYPE"],d["Consequence"],d["MANE"],d["CANONICAL"],d["Amino_acids"],d["Protein_position"])
