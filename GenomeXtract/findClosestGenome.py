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
    Entrez.email = email
    h = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    rec = Entrez.read(h)[0]
    h.close()

    lineage = rec["LineageEx"] + [{"ScientificName": rec["ScientificName"],
                                   "TaxId": rec["TaxId"]}]
    return [(x["ScientificName"], x["TaxId"]) for x in lineage]


# ---------------------------------------------------------
# Nuclear genomes (datasets) — UPDATED WITH OPTION A
# ---------------------------------------------------------

def run_datasets_summary_nuclear(taxid, include_lineage=False):
    """Run datasets summary genome for nuclear genomes, optionally including lineage."""


    cmd = (
        f'datasets summary genome taxon {taxid} '
        f'--assembly-version latest --exclude-atypical --as-json-lines'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        logging.warning(f"datasets command failed:\n{result.stderr}")
        return []

    # Return list of decoded JSON objects
    json_records = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        try:
            json_records.append(json.loads(line))
        except json.JSONDecodeError:
            logging.warning(f"Skipping invalid JSON line:\n{line}")

    #print(json_records)
    return json_records


def parse_nuclear_json(json_records, annotated=False, assembly_level=None):
    """
    Convert datasets JSONL records to a dict:
    { species_name: [ {metadata_dict}, {metadata_dict}, ... ] }
    """
    species_accessions = {}

    for js in json_records:

        # Required basic fields
        accession = js.get("accession")
        org_name = js.get("organism", {}).get("organism_name")

        if not accession or not org_name:
            continue

        asm = js.get("assembly_info", {})
        ann = js.get("annotation_info", {})

        # Filters
        asm_level = asm.get("assembly_level", "").lower()
        if assembly_level and asm_level != assembly_level.lower():
            continue

        annotation_available = bool(ann)
        if annotated and not annotation_available:
            continue

        # Extract meaningful metadata
        metadata = {
            "species": org_name,
            "accession": accession,
            "assembly_name": asm.get("assembly_name", ""),
            "assembly_level": asm.get("assembly_level", ""),
            "assembly_method": asm.get("assembly_method", ""),
            "assembly_date": asm.get("submission_date", ""),
            "bioproject": asm.get("bioproject_accession", ""),
            "biosample": asm.get("biosample_accession", ""),
            "representation": asm.get("representation", ""),
            "chromosomes": asm.get("chromosomes", ""),
            "n50": asm.get("contig_n50", ""),
            "annotation_available": "yes" if annotation_available else "no",
            "provider": ann.get("provider", "") if annotation_available else "",
        }

        species_accessions.setdefault(org_name, []).append(metadata)

    return species_accessions


def list_taxon_nuclear_genomes(taxid, rank, annotated=False, assembly_level=None):
    """
    Option A logic:
    - If rank == species → no lineage
    - If rank != species → include lineage
    """
    include_lineage = (rank != "species")

    json_records = run_datasets_summary_nuclear(
        taxid,
        include_lineage=include_lineage
    )
    return parse_nuclear_json(json_records, annotated, assembly_level)


# ---------------------------------------------------------
# Organellar genomes (plastid / mitochondrial)
# ---------------------------------------------------------

def build_entrez_term_taxid(taxid, genome_type):
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
        return {}

    h2 = Entrez.efetch(db="nucleotide", id=ids, retmode="xml")
    recs = Entrez.read(h2)
    h2.close()

    species_accessions = {}
    for rec in recs:
        org = rec.get("GBSeq_organism")
        acc = rec.get("GBSeq_primary-accession")
        if org and acc:
            species_accessions.setdefault(org, []).append(acc)

    return species_accessions


# ---------------------------------------------------------
# Nearest species search (species input only)
# ---------------------------------------------------------

def nearest_species_with_genome(taxon_name, genome_type, email, annotated=False, assembly_level=None):
    rank, taxid = get_taxid_and_rank(taxon_name, email)
    if not taxid:
        return None

    lineage = get_lineage(taxid, email)

    for name, _tid in reversed(lineage):
        if genome_type in ["chloroplast", "mitochondrial"]:
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
                    acc = rec.get("GBSeq_primary-accession")
                    if acc:
                        return {rec["GBSeq_organism"]: [acc]}

        else:  # nuclear → must use datasets
            r, t = get_taxid_and_rank(name, email)
            sp_dict = list_taxon_nuclear_genomes(t, r, annotated, assembly_level)
            if sp_dict:
                return sp_dict

    return None


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    parser = ArgumentParser(description="Find the closest available reference genome(s) of a given taxon in NCBI.")
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
    # Sanitize taxon name for filesystem (replace spaces with underscores)
    taxon_clean = args.taxon.replace(" ", "_")

    # Build output file path
    out = os.path.join(args.outfolder, f"{taxon_clean}_{args.genome_type}.tsv")


    # ------------------------
    # Species-level search
    # ------------------------
    if rank == "species":
        logging.info("Species rank detected → using nearest-species search.")
        species_dict = nearest_species_with_genome(
            args.taxon,
            args.genome_type,
            args.email,
            args.annotated,
            args.assembly_level
        )
        if species_dict:
            with open(out, "w") as f:
                headers = [
                    "Species", "Accession", "Assembly_Name", "Assembly_Level",
                    "Assembly_Method", "Assembly_Date", "Provider",
                    "Bioproject", "Biosample", "Annotated",
                    "Representation", "N50", "Chromosomes"
                ]
                f.write("\t".join(headers) + "\n")

                for sp, entries in species_dict.items():
                    for md in entries:
                        row = [
                            md.get("species", ""),
                            md.get("accession", ""),
                            md.get("assembly_name", ""),
                            md.get("assembly_level", ""),
                            md.get("assembly_method", ""),
                            md.get("assembly_date", ""),
                            md.get("provider", ""),
                            md.get("bioproject", ""),
                            md.get("biosample", ""),
                            md.get("annotation_available", ""),
                            md.get("representation", ""),
                            str(md.get("n50", "")),
                            str(md.get("chromosomes", ""))
                        ]
                        f.write("\t".join(row) + "\n")

            logging.info(f"Nearest species with a {args.genome_type} saved to {out}")
        else:
            logging.warning("No genome found in the lineage.")
        return
    # ------------------------
    # Higher-level taxon search
    # ------------------------
    logging.info("Higher taxon detected → listing all descendant species with genomes.")

    if args.genome_type in ["chloroplast", "mitochondrial"]:
        species_dict = list_entire_taxon_genomes(taxid, args.genome_type, args.email)
    else:
        species_dict = list_taxon_nuclear_genomes(
            taxid,
            rank,
            args.annotated,
            args.assembly_level
        )

    if species_dict:
        with open(out, "w") as f:
            headers = [
                "Species", "Accession", "Assembly_Name", "Assembly_Level",
                "Assembly_Method", "Assembly_Date", "Provider",
                "Bioproject", "Biosample", "Annotated",
                "Representation", "N50", "Chromosomes"
            ]
            f.write("\t".join(headers) + "\n")

            for sp, entries in species_dict.items():
                for md in entries:
                    row = [
                        md.get("species", ""),
                        md.get("accession", ""),
                        md.get("assembly_name", ""),
                        md.get("assembly_level", ""),
                        md.get("assembly_method", ""),
                        md.get("assembly_date", ""),
                        md.get("provider", ""),
                        md.get("bioproject", ""),
                        md.get("biosample", ""),
                        md.get("annotation_available", ""),
                        md.get("representation", ""),
                        str(md.get("n50", "")),
                        str(md.get("chromosomes", ""))
                    ]
                    f.write("\t".join(row) + "\n")

        logging.info(f"Found {len(species_dict)} species.")
        logging.info(f"Saved list to {out}.")
    else:
        logging.warning("No genomes found for the given taxon.")


if __name__ == "__main__":
    main()