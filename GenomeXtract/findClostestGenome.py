#!/usr/bin/env python3

"""
Find nearest or all available plastome, mitogenome, or nuclear genomes
for a given species or higher taxon using NCBI databases.

This version ensures nuclear genome retrieval works reliably for both
species- and genus-level TaxIDs using NCBI Datasets CLI.

License: GPLv3
"""

import os
import subprocess
import json
import logging
from argparse import ArgumentParser
from Bio import Entrez

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")


# ----------------------------
# Taxonomy helpers
# ----------------------------
def get_taxid_and_rank(taxon_name, email):
    Entrez.email = email
    h = Entrez.esearch(db="taxonomy", term=taxon_name)
    rec = Entrez.read(h)
    h.close()
    if not rec["IdList"]:
        return None, None
    taxid = rec["IdList"][0]

    h2 = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    r = Entrez.read(h2)[0]
    h2.close()
    return r["Rank"], taxid


# ----------------------------
# Nuclear genomes via Datasets CLI
# ----------------------------
def list_taxon_nuclear_genomes(taxid, annotated=False, assembly_level=None):
    """
    List all nuclear genomes for a given TaxID (species or genus).
    Returns list of dicts with species, taxid, accession, genome_type, assembly_level, size.
    """
    level = assembly_level if assembly_level else "all"
    datasets_path = "/Users/kevin/miniconda3/envs/ncbi_datasets/bin/datasets"

    cmd = (
        f'{datasets_path} summary genome taxon {taxid} '
        f'--assembly-level "{level}" --assembly-version latest '
        f'--exclude-atypical --as-json-lines'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        logging.warning(f"Datasets CLI failed for taxid {taxid} or returned no output.")
        return []

    genomes = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        try:
            js = json.loads(line)
            if annotated and js.get("annotation") != "ANNOTATED":
                continue
            genome_info = {
                "species": js.get("organism", {}).get("organism_name"),
                "taxid": js.get("organism", {}).get("taxid"),
                "accession": js.get("assembly", {}).get("accession"),
                "genome_type": "nuclear",
                "assembly_level": js.get("assembly", {}).get("assembly_level"),
                "size": js.get("assembly", {}).get("size", "")
            }
            if genome_info["species"] and genome_info["accession"]:
                genomes.append(genome_info)
        except json.JSONDecodeError:
            continue
    return genomes


# ----------------------------
# Organellar genomes via Entrez
# ----------------------------
def build_entrez_term_taxid(taxid, genome_type):
    term = f"(txid{taxid}[Organism:exp]) AND \"complete genome\""
    if genome_type == "chloroplast":
        term += " AND chloroplast"
    elif genome_type == "mitochondrial":
        term += " AND mitochondrion[Title]"
    return term


def list_entire_taxon_genomes(taxid, genome_type, email):
    Entrez.email = email
    term = build_entrez_term_taxid(taxid, genome_type)
    h = Entrez.esearch(db="nucleotide", term=term, retmax=5000)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        return []

    h2 = Entrez.efetch(db="nucleotide", id=ids, retmode="xml")
    recs = Entrez.read(h2)
    h2.close()

    genomes = []
    for rec in recs:
        genome_info = {
            "species": rec.get("GBSeq_organism"),
            "taxid": rec.get("GBSeq_taxonomy", ""),
            "accession": rec.get("GBSeq_accession-version"),
            "genome_type": genome_type,
            "assembly_level": "complete",
            "size": rec.get("GBSeq_length", "")
        }
        if genome_info["species"] and genome_info["accession"]:
            genomes.append(genome_info)
    return genomes


# ----------------------------
# Main
# ----------------------------
def main():
    parser = ArgumentParser(description="Find nearest or all available genomes in NCBI.")
    parser.add_argument('-g', '--taxon', required=True, help="Species or higher taxon name.")
    parser.add_argument('-o', '--outfolder', required=True, help="Output folder.")
    parser.add_argument('-t', '--genome_type', required=True,
                        choices=['chloroplast', 'mitochondrial', 'nuclear_genome'])
    parser.add_argument('--annotated', action='store_true', help="Only gene-annotated nuclear genomes.")
    parser.add_argument('--assembly_level', required=False,
                        choices=['scaffold', 'chromosome'])
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--email', required=True)
    args = parser.parse_args()

    if os.path.exists(args.outfolder):
        if args.overwrite:
            import shutil
            shutil.rmtree(args.outfolder)
        else:
            raise FileExistsError(f"Output folder {args.outfolder} exists.")
    os.makedirs(args.outfolder, exist_ok=True)

    rank, taxid = get_taxid_and_rank(args.taxon, args.email)
    if not taxid:
        logging.error(f"Taxon {args.taxon} not found in NCBI.")
        return

    logging.info(f"Taxon '{args.taxon}' is rank '{rank}', taxid={taxid}")

    if args.genome_type == "nuclear_genome":
        genomes = list_taxon_nuclear_genomes(taxid, args.annotated, args.assembly_level)
    else:
        genomes = list_entire_taxon_genomes(taxid, args.genome_type, args.email)

    out_file = os.path.join(args.outfolder, f"available_{args.genome_type}_genomes.tsv")
    header = ["species", "taxid", "accession", "genome_type", "assembly_level", "size"]
    with open(out_file, "w") as f:
        f.write("\t".join(header) + "\n")
        for g in genomes:
            f.write("\t".join(str(g.get(h, "")) for h in header) + "\n")

    logging.info(f"Saved {len(genomes)} {args.genome_type} genomes to {out_file}")


if __name__ == "__main__":
    main()
