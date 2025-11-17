#!/usr/bin/env python

"""
This script finds the closest available plastome, mitogenome, and nuclear genome from a given species in the public NCBI database. 

License:
    Copyright 2025 Kevin Karbstein
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import shutil
import subprocess
import json
import logging
from argparse import ArgumentParser
from Bio import Entrez

# Logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")


# ---------------------------------------------------------
# Setup / Validation
# ---------------------------------------------------------

def validate_inputs(taxon, genome_type):
    if not taxon or not genome_type:
        raise ValueError("Both taxon and genome_type must be provided.")


def setup_folder(outfolder, overwrite):
    if os.path.exists(outfolder):
        if overwrite:
            logging.info(f"Overwriting existing folder: {outfolder}")
            shutil.rmtree(outfolder)
        else:
            raise FileExistsError(f"Output folder {outfolder} already exists. Use --overwrite.")
    os.makedirs(outfolder, exist_ok=True)


# ---------------------------------------------------------
# Taxonomy helpers
# ---------------------------------------------------------

def get_taxid_and_rank(taxon_name, email):
    """Return (rank, taxid) for a taxon name."""
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


def get_lineage(taxid, email):
    """Return lineage including self (name, taxid)."""
    Entrez.email = email
    h = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    rec = Entrez.read(h)[0]
    h.close()

    lineage = rec["LineageEx"] + [{"ScientificName": rec["ScientificName"],
                                   "TaxId": rec["TaxId"]}]
    return [(x["ScientificName"], x["TaxId"]) for x in lineage]


# ---------------------------------------------------------
# Nuclear genomes via NCBI Datasets
# ---------------------------------------------------------

def nuclear_genomes_exist(taxon, annotated=False, assembly_level=None):
    """Return True if nuclear genome assemblies exist for taxon."""
    level = assembly_level if assembly_level else "all"

    cmd = (
        f'source ~/.zshrc && conda activate ncbi_datasets && '
        f'datasets summary genome taxon "{taxon}" '
        f'--assembly-level "{level}" --assembly-version latest '
        f'--exclude-atypical --as-json-lines && conda deactivate'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        return False

    assemblies = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        try:
            assemblies.append(json.loads(line))
        except:
            continue

    if annotated:
        assemblies = [a for a in assemblies if a.get("annotation") == "ANNOTATED"]

    return len(assemblies) > 0


def list_taxon_nuclear_genomes(taxid, annotated=False, assembly_level=None):
    """List all nuclear genomes for ANY descendant of a taxid."""
    level = assembly_level if assembly_level else "all"

    cmd = (
        f'source ~/.zshrc && conda activate ncbi_datasets && '
        f'datasets summary genome taxon {taxid} '
        f'--assembly-level "{level}" --assembly-version latest '
        f'--exclude-atypical --as-json-lines && conda deactivate'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return []

    species = set()

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        try:
            js = json.loads(line)
            org = js.get("organism", {}).get("organism_name", None)
            if not org:
                continue
            if annotated and js.get("annotation") != "ANNOTATED":
                continue
            species.add(org)
        except:
            continue

    return sorted(species)


# ---------------------------------------------------------
# Organellar genomes (plastid / mitochondrial)
# ---------------------------------------------------------

def build_entrez_term_taxid(taxid, genome_type):
    """Taxon-expanded search using txidXXXX[Organism:exp]."""
    term = f"(txid{taxid}[Organism:exp]) AND \"complete genome\""

    if genome_type == "chloroplast":
        term += " AND chloroplast"
    elif genome_type == "mitochondrial":
        term += " AND mitochondrion[Title]"
    elif genome_type == "nuclear_genome":
        term += " AND genomic DNA[filter]"
    else:
        raise ValueError(f"Invalid genome_type: {genome_type}")

    return term


def list_entire_taxon_genomes(taxid, genome_type, email):
    """List all species with an organellar genome inside descendant taxa."""
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

    species = set()
    for rec in recs:
        org = rec.get("GBSeq_organism")
        if org:
            species.add(org)

    return sorted(species)


# ---------------------------------------------------------
# Nearest-species logic (species input only)
# ---------------------------------------------------------

def nearest_species_with_genome(taxon_name, genome_type, email, annotated=False, assembly_level=None):
    rank, taxid = get_taxid_and_rank(taxon_name, email)
    if not taxid:
        return None

    lineage = get_lineage(taxid, email)

    for name, _tid in reversed(lineage):

        if genome_type in ["chloroplast", "mitochondrial"]:
            # Organellar genomes via Entrez
            term = f"\"{name}\"[Organism] AND \"complete genome\""
            if genome_type == "chloroplast":
                term += " AND chloroplast"
            else:
                term += " AND mitochondrion[Title]"

            h = Entrez.esearch(db="nucleotide", term=term, retmax=20)
            ids = Entrez.read(h)["IdList"]
            h.close()

            if ids:
                h2 = Entrez.efetch(db="nucleotide", id=ids, retmode="xml")
                recs = Entrez.read(h2)
                h2.close()

                for rec in recs:
                    return rec["GBSeq_organism"]

        else:  # nuclear
            if nuclear_genomes_exist(name, annotated, assembly_level):
                return name

    return None


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    parser = ArgumentParser(description="Find the closest available reference genomes(s) of a given taxon in NCBI.")
    parser.add_argument('-g', '--taxon', required=True, help="Species or higher-level taxon name (e.g., Genus or Family).")
    parser.add_argument('-o', '--outfolder', required=True, help="Output folder for the result file.")
    parser.add_argument('-t', '--genome_type', required=True,
                        choices=['chloroplast', 'mitochondrial', 'nuclear_genome'])
    parser.add_argument('--annotated', action='store_true',
                        help="Select only gene-annotated nuclear genomes.")
    parser.add_argument('--assembly_level', required=False,
                        choices=['scaffold', 'chromosome'],
                        help="Choose the assembly level of the nuclear genome (STRING; scaffold, chromosome).")
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--email', required=True)
    args = parser.parse_args()

    validate_inputs(args.taxon, args.genome_type)
    setup_folder(args.outfolder, args.overwrite)

    rank, taxid = get_taxid_and_rank(args.taxon, args.email)

    if not taxid:
        logging.error("Taxon not found in NCBI.")
        return

    logging.info(f"Taxon '{args.taxon}' is rank '{rank}', taxid={taxid}")

    #
    # CASE 1: species → nearest-species behavior
    #
    if rank == "species":
        logging.info("Species rank detected → using nearest-species search.")

        species = nearest_species_with_genome(
            args.taxon,
            args.genome_type,
            args.email,
            args.annotated,
            args.assembly_level
        )

        out = os.path.join(args.outfolder, "result.txt")

        if species:
            logging.info(f"Nearest species with a {args.genome_type}: {species}")
            with open(out, "w") as f:
                f.write(species + "\n")
        else:
            logging.warning("No genome found in the lineage.")
        return

    #
    # CASE 2: higher taxon → list all descendant species
    #
    logging.info("Higher taxon detected → listing all descendant species with genomes.")

    if args.genome_type in ["chloroplast", "mitochondrial"]:
        species_list = list_entire_taxon_genomes(
            taxid, args.genome_type, args.email
        )
    else:
        species_list = list_taxon_nuclear_genomes(
            taxid, args.annotated, args.assembly_level
        )

    out = os.path.join(args.outfolder, "available_species.txt")
    with open(out, "w") as f:
        f.write("\n".join(species_list))

    logging.info(f"Found {len(species_list)} species.")
    logging.info(f"Saved list to {out}.")


if __name__ == "__main__":
    main()