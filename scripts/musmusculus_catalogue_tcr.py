from pathlib import Path
import json

alleles = dict()

with open(Path("data") / "musmusculus_tcr.fasta", "r") as f:
    for line in f.readlines():
        if line.startswith(">"):
            fields = line.split("|")
            allele_name = fields[1]
            gene = allele_name.split("*")[0]
            allele_designation = allele_name.split("*")[1]
            functionality = fields[3].strip("()[]")

            if not gene in alleles:
                alleles[gene] = dict()

            alleles[gene][allele_designation] = functionality


with open(
    Path("src") / "tidytcells" / "_resources" / "valid_musmusculus_tr.json", "w"
) as f:
    json.dump(alleles, f, indent=4)
