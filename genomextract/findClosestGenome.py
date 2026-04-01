#!/usr/bin/env python

"""
This script searches the public NCBI database for the closest available plastid, mitochondrial or nuclear genome sequence(s) from a given species or taxonomic rank. 
It ranks the species according to their genetic similarity to the target taxon sequence based on alignment average nucleotide identity (ANI) for 
plastid genomes or by using mash-based distance for mitochondrial and nuclear genomes. Alternatively, the script can be used to find all genome sequences 
for a given taxonomic group by selecting the largest genome sequence or a user-specfied genome size as target. The script also automatically filters misaligned samples (organellar genomes), 
or collapses all individuals of a given taxon for genetic similarity comparisons. 


License:
    Copyright 2026 Kevin Karbstein
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
from Bio import SeqIO
from Bio.Entrez.Parser import StringElement
import tempfile
import time
import numpy as np
import glob
import pandas as pd
from Bio.Seq import Seq



# -------------------
# Setup / Validation
# -------------------

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


def get_tmpdir(outfolder):
    """
    Create a temporary working directory inside the output folder.
    """
    tmp_base = os.path.join(outfolder, "tmp")
    os.makedirs(tmp_base, exist_ok=True)

    return tempfile.mkdtemp(dir=tmp_base)


# -------------------------------------------
# Find and download nuclear genome sequences
# -------------------------------------------

#def find_downloaded_genome_fastas(base_dir, accessions, acc2species):
#    genome_files = {}
#
#    for root, _, files in os.walk(base_dir):
#        for fname in files:
#            if not fname.endswith(".fna"):
#                continue
#
#            full_path = os.path.join(root, fname)
#
#            matched_acc = None
#            for acc in accessions:
#                if acc in full_path:
#                    matched_acc = acc
#                    break
#
#            if matched_acc is None:
#                continue
#
#            if matched_acc not in genome_files:
#                genome_files[matched_acc] = {
#                    "path": full_path,
#                    "species": acc2species.get(matched_acc, "unknown")
#                }
#
#    return genome_files


def download_nuclear_genomes_batch(
    accessions,
    species_dict,
    outdir,
    dehydrated=True,
    chunk_size=6,
    retries=4,
    include="genome",
):
    """
    Download nuclear genomes in chunks via NCBI datasets.

    Returns:
        {
            accession: {
                "path": "/path/to/file.fna",
                "species": "Species name"
            }
        }

    Behavior:
    - tolerates partial failures during rehydrate
    - continues with successfully downloaded FASTA files
    - logs which accessions are missing after each chunk
    """
    if not accessions:
        logging.warning("No accessions provided for download.")
        return {}

    accessions = sorted(set(accessions))
    acc2species = build_accession_to_species_map(species_dict)
    genome_files = {}

    all_chunks = list(chunk_list(accessions, chunk_size))
    logging.info(
        f"Downloading {len(accessions)} nuclear genomes in {len(all_chunks)} chunk(s)"
    )

    for chunk_idx, chunk in enumerate(all_chunks, start=1):
        chunk_label = f"chunk {chunk_idx}/{len(all_chunks)}"
        chunk_dir = os.path.join(outdir, f"batch_{chunk_idx:03d}")
        os.makedirs(chunk_dir, exist_ok=True)

        acc_file = os.path.join(chunk_dir, "accessions.txt")
        out_zip = os.path.join(chunk_dir, "nuclear_genomes.zip")
        extract_dir = os.path.join(chunk_dir, "nuclear_genomes")

        with open(acc_file, "w") as f:
            for acc in chunk:
                f.write(acc + "\n")

        cmd = (
            f'datasets download genome accession '
            f'--inputfile "{acc_file}" '
            f'--filename "{out_zip}" '
            f'--include {include} '
            f'--no-progressbar'
        )

        if dehydrated:
            cmd += " --dehydrated"

        # -------------------------
        # Step 1: download zip
        # -------------------------
        try:
            run_cmd_with_retries(
                cmd,
                retries=retries,
                wait_seconds=15,
                label=f"datasets download {chunk_label}"
            )
        except RuntimeError as e:
            logging.warning(f"Skipping {chunk_label}: datasets download failed\n{e}")
            continue

        # -------------------------
        # Step 2: unzip
        # -------------------------
        os.makedirs(extract_dir, exist_ok=True)
        try:
            run_cmd_with_retries(
                f'unzip -oq "{out_zip}" -d "{extract_dir}"',
                retries=2,
                wait_seconds=3,
                label=f"unzip {chunk_label}"
            )
        except RuntimeError as e:
            logging.warning(f"Skipping {chunk_label}: unzip failed\n{e}")
            continue

        # -------------------------
        # Step 3: rehydrate
        # -------------------------
        if dehydrated:
            try:
                run_cmd_with_retries(
                    f'datasets rehydrate --directory "{extract_dir}"',
                    retries=retries,
                    wait_seconds=20,
                    label=f"rehydrate {chunk_label}"
                )
            except RuntimeError as e:
                logging.warning(
                    f"Partial rehydrate failure in {chunk_label}; continuing with available files.\n{e}"
                )

        # -------------------------
        # Step 4: collect available FASTA files
        # -------------------------
        chunk_files = find_downloaded_genome_fastas(extract_dir, chunk, acc2species)

        if not chunk_files:
            logging.warning(f"No .fna files found after processing {chunk_label}")
            continue

        genome_files.update(chunk_files)

        downloaded_accs = set(chunk_files.keys())
        missing_accs = set(chunk) - downloaded_accs

        logging.info(
            f"{chunk_label}: recovered FASTA for {len(downloaded_accs)}/{len(chunk)} accession(s)"
        )

        if missing_accs:
            logging.warning(
                f"{chunk_label}: missing FASTA for {len(missing_accs)} accession(s): "
                + ", ".join(sorted(missing_accs))
            )

    logging.info(
        f"Final nuclear genome download summary: recovered {len(genome_files)}/{len(accessions)} accession(s)"
    )

    if not genome_files:
        logging.error("No nuclear genome FASTA files were recovered.")
    else:
        recovered = set(genome_files.keys())
        missing_total = set(accessions) - recovered
        if missing_total:
            logging.warning(
                f"Overall missing FASTA for {len(missing_total)} accession(s): "
                + ", ".join(sorted(missing_total))
            )

    return genome_files



def find_downloaded_genome_fastas(base_dir, accessions, acc2species):
    """
    Find downloaded genomic FASTA files (*.fna) for requested accessions.

    Returns:
        {
            accession: {
                "path": "/path/to/file.fna",
                "species": "Species name"
            }
        }
    """
    genome_files = {}
    accession_set = set(accessions)

    for root, _, files in os.walk(base_dir):
        for fname in files:
            if not fname.endswith(".fna"):
                continue

            matched_acc = None
            for acc in accession_set:
                if acc in fname or acc in root:
                    matched_acc = acc
                    break

            if matched_acc is None:
                continue

            if matched_acc not in genome_files:
                genome_files[matched_acc] = {
                    "path": os.path.join(root, fname),
                    "species": acc2species.get(matched_acc, "unknown")
                }

    return genome_files


def run_datasets_summary_nuclear(taxid, include_lineage=False, retries=3, retry_wait=5):
    cmd = (
        f"datasets summary genome taxon {taxid} "
        f"--assembly-version latest --exclude-atypical --as-json-lines"
    )

    last_stderr = ""

    for attempt in range(1, retries + 1):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            json_records = []
            for line in result.stdout.strip().split("\n"):
                if not line.strip():
                    continue
                try:
                    json_records.append(json.loads(line))
                except json.JSONDecodeError:
                    logging.warning(f"Skipping invalid JSON line:\n{line}")
            return json_records

        last_stderr = (result.stderr or "").strip()

        #if "no such host" in last_stderr.lower():
            #logging.warning(
                #f"datasets attempt {attempt}/{retries} failed due to DNS/network problem:\n{last_stderr}"
            #)
        #else:
            #logging.warning(f"datasets attempt {attempt}/{retries} failed:\n{last_stderr}")

        if attempt < retries:
            time.sleep(retry_wait)

    raise RuntimeError(
        f"datasets summary genome failed after {retries} attempts.\n{last_stderr}"
    )


def parse_nuclear_json(json_records, annotated=False, assembly_level=None):
    """
    Convert datasets JSONL records to:
    { species_name: [ {metadata_dict}, {metadata_dict}, ... ] }
    """
    species_accessions = {}

    for js in json_records:
        accession = js.get("accession")
        org_name = js.get("organism", {}).get("organism_name")

        if not accession or not org_name:
            continue

        asm = js.get("assembly_info", {})
        stats = js.get("assembly_stats", {})
        ann = js.get("annotation_info", {})

        asm_level_value = asm.get("assembly_level", "")
        if assembly_level and asm_level_value.lower() != assembly_level.lower():
            continue

        annotation_available = bool(ann)
        if annotated and not annotation_available:
            continue

        metadata = {
            "species": org_name,
            "accession": accession,
            "assembly_name": asm.get("assembly_name", ""),
            "assembly_level": asm.get("assembly_level", ""),
            "assembly_method": asm.get("assembly_method", ""),
            "assembly_date": asm.get("release_date", ""),
            "bioproject": asm.get("bioproject_accession", ""),
            "biosample": asm.get("biosample", {}).get("accession", ""),
            "chromosomes": stats.get("total_number_of_chromosomes", ""),
            "n50": stats.get("contig_n50", ""),
            "genome_size": stats.get("total_sequence_length", ""),
            "annotation_available": "yes" if annotation_available else "no",
            "provider": ann.get("provider", "") if annotation_available else "",
        }

        species_accessions.setdefault(org_name, []).append(metadata)

    return species_accessions


def list_taxon_nuclear_genomes(taxid, rank, annotated=False, assembly_level=None):
    """
    Logic:
    - If rank == species → no lineage
    - If rank != species → include lineage
    """
    json_records = run_datasets_summary_nuclear(taxid)
    return parse_nuclear_json(json_records, annotated=annotated, assembly_level=assembly_level)


# ------------------------------------------------------------------------
# Find and download organellar genome sequences (plastid / mitochondrial)
# ------------------------------------------------------------------------

def fetch_entrez_records(ids, batch_size=200):
    all_records = []

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]

        handle = Entrez.efetch(
            db="nucleotide",
            id=batch,
            retmode="xml"
        )
        recs = Entrez.read(handle)
        handle.close()

        all_records.extend(recs)

    return all_records


def fetch_with_retry(ids, retries=3, wait_seconds=3, batch_size=200):
    """
    Fetch Entrez nucleotide records with retries.
    Logs transient issues gently and reports success after retry.
    """
    last_err = None

    for attempt in range(1, retries + 1):
        try:
            records = fetch_entrez_records(ids, batch_size=batch_size)

            if attempt > 1:
                logging.info(
                    f"Entrez efetch succeeded on retry {attempt}/{retries} "
                    f"for {len(ids)} id(s)"
                )

            return records

        except Exception as e:
            last_err = e

            if attempt < retries:
                log_fn = logging.info if attempt == 1 else logging.warning
                log_fn(
                    f"Transient Entrez efetch issue ({attempt}/{retries})\n"
                    f"n_ids={len(ids)}\n"
                    f"error={e}"
                )
                time.sleep(wait_seconds * attempt)
            else:
                logging.error(
                    f"Entrez efetch failed after {retries} attempts\n"
                    f"n_ids={len(ids)}\n"
                    f"error={e}"
                )

    raise RuntimeError(f"Entrez fetch failed after {retries} attempts: {last_err}")


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
    Entrez.email = email
    term = build_entrez_term_taxid(taxid, genome_type)

    rec = entrez_search_with_retry(
        db="nucleotide",
        term=term,
        retmax=500,
        retries=5,
        wait_seconds=3
    )
    ids = rec.get("IdList", [])

    if not ids:
        return {}

    recs = fetch_with_retry(ids)

    species_accessions = {}
    for rec in recs:
        org = str(rec.get("GBSeq_organism", "")).strip()
        acc = str(rec.get("GBSeq_primary-accession", "")).strip()

        acc_taxid = None
        for feature in rec.get("GBSeq_feature-table", []):
            for qual in feature.get("GBFeature_quals", []):
                if qual.get("GBQualifier_name") == "db_xref":
                    val = qual.get("GBQualifier_value", "")
                    if val.startswith("taxon:"):
                        acc_taxid = val.replace("taxon:", "")
                        break
            if acc_taxid:
                break

        if not acc_taxid:
            try:
                _, acc_taxid = get_taxid_and_rank(org, email)
            except Exception:
                acc_taxid = ""

        references = rec.get("GBSeq_references", [])
        provider = ""
        if references:
            first_ref = references[0]
            authors = first_ref.get("GBReference_authors", [])
            if authors:
                provider = ", ".join(authors)

        length = str(rec.get("GBSeq_length", "")).strip()

        if org and acc:
            species_accessions.setdefault(org, []).append({
                "species": org,
                "accession": acc,
                "assembly_name": rec.get("GBSeq_definition", ""),
                "assembly_level": rec.get("GBSeq_topology", ""),
                "assembly_date": rec.get("GBSeq_update-date", ""),
                "provider": provider,
                "bioproject": rec.get("GBSeq_project", ""),
                "taxonomy": rec.get("GBSeq_taxonomy", ""),
                "sequence_length": length,
                "accession_taxid": acc_taxid
            })

    return species_accessions




def download_fasta_entrez(accessions, outdir, email):
    """
    Download FASTA files from NCBI nucleotide for organellar ANI calculations.
    """
    Entrez.email = email
    fasta_files = {}

    for acc in accessions:
        out_fa = os.path.join(outdir, f"{acc}.fasta")

        if os.path.exists(out_fa):
            fasta_files[acc] = out_fa
            continue

        handle = Entrez.efetch(
            db="nucleotide",
            id=acc,
            rettype="fasta",
            retmode="text"
        )

        with open(out_fa, "w") as f:
            f.write(handle.read())

        handle.close()
        fasta_files[acc] = out_fa

    return fasta_files



def chunk_list(items, chunk_size):
    for i in range(0, len(items), chunk_size):
        yield items[i:i + chunk_size]


def run_cmd_with_retries(cmd, retries=4, wait_seconds=10, label="command"):
    last_err = None

    for attempt in range(1, retries + 1):
        logging.info(f"{label}: attempt {attempt}/{retries}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            return result

        last_err = (result.stderr or result.stdout or "").strip()
        logging.warning(f"{label} failed on attempt {attempt}:\n{last_err}")

        if attempt < retries:
            time.sleep(wait_seconds * attempt)

    raise RuntimeError(f"{label} failed after {retries} attempts.\n{last_err}")



def entrez_search_with_retry(
    db,
    term,
    retmax=500,
    retries=5,
    wait_seconds=3,
    suppress_transient_logs=True
):
    last_err = None

    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=retmax)
            result = Entrez.read(handle)
            handle.close()
            return result

        except Exception as e:
            last_err = e

            if attempt < retries:
                if not suppress_transient_logs:
                    logging.warning(
                        f"Transient Entrez esearch failure ({attempt}/{retries})\n"
                        f"db={db}\n"
                        f"term={term}\n"
                        f"error={e}"
                    )
                time.sleep(wait_seconds * attempt)
            #else:
                #logging.error(
                    #f"Entrez esearch failed after {retries} attempts\n"
                    #f"db={db}\n"
                    #f"term={term}\n"
                    #f"error={e}"
                #)

    raise RuntimeError(f"Entrez esearch failed after {retries} attempts: {last_err}")




from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.error import HTTPError
from http.client import IncompleteRead

def download_fasta_parallel(accessions, outdir, email, max_workers=3, retries=4, wait_seconds=3):
    """
    Parallel download FASTA files for organellar genomes.
    Keep concurrency low to avoid NCBI rate limiting.
    """
    Entrez.email = email

    # hard cap for NCBI
    if max_workers is None:
        max_workers = 3
    max_workers = min(max_workers, 3)

    def worker(acc):
        for attempt in range(1, retries + 1):
            try:
                result = download_fasta_entrez([acc], outdir, email)
                return acc, result.get(acc, None)

            except HTTPError as e:
                if e.code == 429:
                    logging.warning(
                        f"{acc} failed ({attempt}/{retries}): HTTP 429 Too Many Requests"
                    )
                else:
                    logging.warning(
                        f"{acc} failed ({attempt}/{retries}): HTTP Error {e.code}"
                    )

            except IncompleteRead as e:
                logging.warning(
                    f"{acc} failed ({attempt}/{retries}): IncompleteRead({len(e.partial)} bytes read)"
                )

            except Exception as e:
                logging.warning(f"{acc} failed ({attempt}/{retries}): {e}")

            if attempt < retries:
                time.sleep(wait_seconds * attempt)

        return acc, None

    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(worker, acc) for acc in accessions]

        for future in as_completed(futures):
            acc, path = future.result()
            if path:
                results[acc] = path

    logging.info(f"Successfully downloaded {len(results)}/{len(accessions)} FASTA files")
    return results


msh_cache = {}
genome_size_cache = {}



# ------------------
# Taxonomy helpers 
# ------------------
def get_taxid_and_rank(taxon_name, email, retries=4, wait_seconds=3):
    """
    Get NCBI TaxID and rank for a given taxon name.

    Retries on:
    - Entrez exceptions
    - empty IdList results (transient backend issues)

    Returns:
        (rank, taxid)
    """
    Entrez.email = email
    taxon_name = taxon_name.strip()

    record = None
    last_err = None

    # -------------------------
    # Step 1: esearch with retries
    # -------------------------
    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.esearch(
                db="taxonomy",
                term=f'"{taxon_name}"[Scientific Name]'
            )
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])

            if id_list:
                break

            # empty IdList can be transient -> retry
            logging.warning(
                f"Entrez taxonomy search returned no taxid ({attempt}/{retries}) "
                f"for '{taxon_name}'"
            )
            record = None

        except Exception as e:
            last_err = e
            logging.warning(
                f"Entrez taxonomy search returned empty result ({attempt}/{retries}) "
                f"for '{taxon_name}' - retrying"
            )
            record = None

        if attempt < retries:
            time.sleep(wait_seconds * attempt)

    if record is None or not record.get("IdList"):
        if last_err is not None:
            raise RuntimeError(f"Entrez search failed for '{taxon_name}': {last_err}")
        raise ValueError(f"No taxid found for '{taxon_name}' after {retries} attempts")

    taxid = record["IdList"][0]

    # -------------------------
    # Step 2: efetch with retries
    # -------------------------
    rec = None
    last_err = None

    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            fetched = Entrez.read(handle)
            handle.close()

            if fetched and len(fetched) > 0:
                rec = fetched[0]
                break

            logging.warning(
                f"Entrez taxonomy fetch returned empty result ({attempt}/{retries}) "
                f"for taxid {taxid}"
            )

        except Exception as e:
            last_err = e
            logging.warning(
                f"Entrez taxonomy fetch failed ({attempt}/{retries}) for "
                f"taxid {taxid}: {e}"
            )

        if attempt < retries:
            time.sleep(wait_seconds * attempt)

    if rec is None:
        if last_err is not None:
            raise RuntimeError(f"Entrez fetch failed for taxid {taxid}: {last_err}")
        raise RuntimeError(f"Entrez fetch returned no data for taxid {taxid}")

    rank = rec.get("Rank", "no_rank")
    taxid = rec.get("TaxId", taxid)

    return rank, taxid


def select_largest_entry(entries, genome_type):
    """
    Select the largest genome entry from a list of metadata dicts.
    For organellar genomes uses 'sequence_length'.
    For nuclear genomes uses genome_size as fallback.
    """
    if not entries:
        return None

    if genome_type in ["chloroplast", "mitochondrial"]:
        return max(entries, key=lambda g: int(g.get("sequence_length", 0) or 0))
    else:
        return max(entries, key=lambda g: int(g.get("genome_size", 0) or 0))


def get_lineage_with_ranks(taxid, email):
    """
    Return full lineage including the queried taxon itself:
    [
        {"name": "...", "taxid": "...", "rank": "..."},
        ...
    ]
    """
    Entrez.email = email
    h = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    rec = Entrez.read(h)[0]
    h.close()

    lineage = rec.get("LineageEx", [])
    lineage.append({
        "ScientificName": rec.get("ScientificName", ""),
        "TaxId": rec.get("TaxId", ""),
        "Rank": rec.get("Rank", "no_rank")
    })

    result = []
    for x in lineage:
        result.append({
            "name": x.get("ScientificName", ""),
            "taxid": str(x.get("TaxId", "")),
            "rank": x.get("Rank", "no_rank")
        })

    return result


def get_search_ranks_upto(max_rank, base_ranks=("species", "genus", "family", "order")):
    """
    Return all ranks from the start of base_ranks up to and including max_rank.
    """
    if max_rank not in base_ranks:
        raise ValueError(f"Invalid rank: {max_rank}")

    idx = base_ranks.index(max_rank)
    return base_ranks[:idx + 1]


def trim_ranks_from_input_rank(rank, ranks):
    """
    Trim a rank list so that searching starts at the user's input rank.
    """
    if rank not in ranks:
        return ranks

    idx = ranks.index(rank)
    return ranks[idx:]


def extract_rank_targets_from_lineage(lineage, ranks):
    """
    Extract selected ranks from lineage in the requested order.

    Returns:
        [
            {"rank": "species", "name": "...", "taxid": "..."},
            {"rank": "genus", "name": "...", "taxid": "..."},
            ...
        ]
    """
    rank_map = {}

    for item in lineage:
        r = item.get("rank")
        if r in ranks:
            rank_map[r] = {
                "rank": r,
                "name": item.get("name"),
                "taxid": item.get("taxid")
            }

    return [rank_map[r] for r in ranks if r in rank_map]



def get_genome_size(path):
    """
    Get genome size with caching to avoid repeated disk reads.
    """
    if path in genome_size_cache:
        return genome_size_cache[path]

    size = sum(len(r.seq) for r in SeqIO.parse(path, "fasta"))
    genome_size_cache[path] = size
    return size


def select_target_genome(
    genome_files,
    taxon_name,
    genome_size_mb=None,
    size_tolerance=0.2
):
    """
    Select target genome accession.

    Strategy:
    1. exact species match
    2. genus match + Mash medoid
    3. size-based fallback + Mash medoid
    4. Mash medoid across all genomes
    """
    target_name = taxon_name.strip().lower()

    if genome_size_mb is not None:
        query_size = genome_size_mb * 1e6
        logging.info(f"Using user-defined genome size: {genome_size_mb} Mb")
    else:
        query_size = None

    # 1. Exact species match
    exact_matches = {
        acc: meta for acc, meta in genome_files.items()
        if meta["species"].strip().lower() == target_name
    }
    if exact_matches:
        if len(exact_matches) == 1:
            best_acc = next(iter(exact_matches))
            logging.info(f"Using target genome (exact species match): {best_acc}")
            return best_acc
        best_acc, dist = mash_select_best(exact_matches)
        logging.info(f"Selected exact-species representative: {best_acc} (dist={dist})")
        return best_acc

    # 2. Genus match
    genus_matches = filter_by_genus(genome_files, taxon_name)
    if genus_matches:
        logging.info(f"Found {len(genus_matches)} genus-level candidates")
        best_acc, dist = mash_select_best(genus_matches)
        logging.info(f"Selected (genus + Mash): {best_acc} (dist={dist})")
        return best_acc

    logging.warning("No genus match found -> trying size-based fallback")

    # 3. Size-based fallback
    if query_size is not None:
        size_matches = filter_by_size(genome_files, query_size, tolerance=size_tolerance)
        if size_matches:
            logging.info(f"Found {len(size_matches)} size-matched candidates")
            best_acc, dist = mash_select_best(size_matches)
            logging.info(f"Selected (size + Mash): {best_acc} (dist={dist})")
            return best_acc

    # 4. Global fallback
    logging.warning("No size match or no size provided -> selecting global Mash medoid")
    best_acc, dist = mash_select_best(genome_files)
    logging.info(f"Selected (global Mash medoid): {best_acc} (dist={dist})")
    return best_acc


def collect_genomes_across_ranks(
    taxon_name,
    taxid,
    genome_type,
    email,
    annotated=False,
    assembly_level=None,
    ranks=("species", "genus"),
    max_genomes=None
):
    lineage = get_lineage_with_ranks(taxid, email)
    rank_targets = extract_rank_targets_from_lineage(lineage, ranks=ranks)

    logging.info(
        "Rank search order: " +
        " -> ".join([f"{x['rank']}:{x['name']}" for x in rank_targets])
    )

    combined = {}
    seen_accessions = set()
    failed_ranks = []

    for target in rank_targets:
        r = target["rank"]
        t = target["taxid"]
        n = target["name"]

        logging.info(f"Searching {genome_type} genomes at rank={r}, taxon={n}, taxid={t}")

        try:
            if genome_type in ["chloroplast", "mitochondrial"]:
                rank_dict = list_entire_taxon_genomes(t, genome_type, email)
            else:
                rank_dict = list_taxon_nuclear_genomes(
                    t,
                    r,
                    annotated=annotated,
                    assembly_level=assembly_level
                )
        except RuntimeError as e:
            logging.warning(
                f"Skipping rank={r}, taxon={n} because query failed:\n{e}"
            )
            failed_ranks.append((r, n))
            continue

        if not rank_dict:
            logging.info(f"No genomes found at rank={r}, taxon={n}")
            continue

        added_here = 0

        for species, entries in rank_dict.items():
            for entry in entries:
                acc = entry.get("accession")
                if not acc or acc in seen_accessions:
                    continue

                combined.setdefault(species, []).append(entry)
                seen_accessions.add(acc)
                added_here += 1

                if max_genomes is not None and len(seen_accessions) >= max_genomes:
                    logging.info(f"Reached max_genomes={max_genomes}")
                    return combined

        logging.info(
            f"Added {added_here} accession(s) from rank={r}, taxon={n}. "
            f"Total unique accessions so far: {len(seen_accessions)}"
        )

    if not combined and failed_ranks:
        raise RuntimeError(
            "No genomes could be collected because all rank queries failed: " +
            ", ".join([f"{r}:{n}" for r, n in failed_ranks])
    )

    if failed_ranks:
        logging.warning(
            "Some rank queries failed but the search continued: " +
            ", ".join([f"{r}:{n}" for r, n in failed_ranks])
        )

    return combined



def filter_by_genus(genome_files, taxon_name):
    """
    Filter genomes by genus using stored species names.
    """
    genus = taxon_name.split()[0].lower()

    return {
        acc: meta
        for acc, meta in genome_files.items()
        if meta["species"].lower().startswith(genus + " ")
    }


def filter_by_size(genome_files, query_size, tolerance=0.5):
    """
    Keep genomes within ± tolerance of query size.
    """
    filtered = {}

    for acc, meta in genome_files.items():
        size = get_genome_size(meta["path"])
        if abs(size - query_size) / query_size <= tolerance:
            filtered[acc] = meta

    return filtered





# -------------------------
# Determine target ID
# -------------------------
def determine_target_id(fasta_files, genome_size_mb=None):
    """
    Determine a stable target_id for ANI calculation.

    Strategy:
    - If genome_size_mb is provided → choose genome closest in size
    - Otherwise → choose longest genome (default fallback)
    """

    genome_sizes = {}

    # -------------------------
    # Calculate genome sizes
    # -------------------------
    for acc, fasta in fasta_files.items():
        size = sum(len(r.seq) for r in SeqIO.parse(fasta, "fasta"))
        genome_sizes[acc] = size
        

    if not genome_sizes:
        raise RuntimeError("Could not determine genome sizes.")

    # -------------------------
    # OPTION 1: User-defined size
    # -------------------------
    if genome_size_mb is not None:
        target_size = genome_size_mb * 1e6

        best_acc = min(
            genome_sizes,
            key=lambda acc: abs(genome_sizes[acc] - target_size)
        )

        logging.info(
            f"[INFO] No exact genome found for target taxon → "
            f"selecting proxy genome based on user-defined size ({genome_size_mb} Mb)"
        )

        logging.info(
            f"[INFO] Selected target_id: {best_acc} "
            f"(size={genome_sizes[best_acc]} bp)"
        )

        return best_acc

    # -------------------------
    # OPTION 2: Default fallback
    # -------------------------
    best_acc = max(genome_sizes, key=genome_sizes.get)

    logging.warning(
        "No genome available for the requested taxon.\n"
        "→ Using proxy genome for ANI calculation.\n"
        "→ Selection strategy: longest available genome.\n"
        "You can override this using --genome-size-mb."
    )

    logging.info(
        f"[INFO] Selected target_id: {best_acc} "
        f"(size={genome_sizes[best_acc]} bp)"
    )

    return best_acc



def build_accession_to_species_map(species_dict):
    acc2species = {}

    for species, entries in species_dict.items():
        for entry in entries:
            acc = entry.get("accession")
            if acc:
                acc2species[acc] = species

    return acc2species



# --------------------------------------
# Create alignment (chloroplast genome)
# --------------------------------------

def run_mafft(fasta_files, output_alignment, taxon_name=None, outfolder=None):
    """
    Create plastid genome alignment to calculate ANI.
    Expects exactly one FASTA record per input file.
    Optionally saves a copy as: taxon_mafft_alignment.fasta
    """

    input_concat = output_alignment + ".input.fasta"

    # -------------------------
    # Create input fasta
    # -------------------------
    with open(input_concat, "w") as out:
        for acc, fasta_path in fasta_files.items():
            records = list(SeqIO.parse(fasta_path, "fasta"))

            if len(records) == 0:
                raise RuntimeError(f"No FASTA record found in file: {fasta_path}")

            if len(records) > 1:
                raise RuntimeError(
                    f"Expected exactly 1 FASTA record in {fasta_path}, "
                    f"but found {len(records)}. "
                    f"Please ensure one genome per accession file."
                )

            rec = records[0]
            rec.id = acc
            rec.name = acc
            rec.description = ""
            SeqIO.write(rec, out, "fasta")

    # -------------------------
    # Execute MAFFT
    # -------------------------
    with open(output_alignment, "w") as out_handle:
        
        n_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))

        subprocess.run(
            ["mafft", "--thread", str(n_threads), "--adjustdirection", "--auto", input_concat],
            check=True,
            stdout=out_handle,
            text=True
        )

    # -------------------------
    # Save the alignment
    # -------------------------
    if taxon_name and outfolder:
        taxon_clean = taxon_name.replace(" ", "_")
        final_path = os.path.join(outfolder, f"{taxon_clean}_mafft_alignment.fasta")
        shutil.copy(output_alignment, final_path)
        logging.info(f"Saved MAFFT alignment to: {final_path}")

    return output_alignment


# -----------------------------------------
# Filter MAFFT Aligmment for missing data
# -----------------------------------------

def filter_alignment_by_shared_coverage(
    alignment_file,
    output_file,
    min_shared_fraction=0.7,
    min_site_occupancy=0.5,
    force_keep_ids=None
):
    """
    Remove sequences with low shared coverage and remove poorly occupied columns.

    Parameters:
        min_shared_fraction:
            Fraction of well-supported columns in which a sequence must have data.
        min_site_occupancy:
            Fraction of sequences that must have a non-gap/non-N base at a site
            for the site to be kept.
    """

    records = list(SeqIO.parse(alignment_file, "fasta"))

    if not records:
        raise ValueError("Alignment is empty")

    aln = np.array([list(str(r.seq).upper()) for r in records])
    n_seq, n_sites = aln.shape

    # -------------------------
    # Step 1: determine good columns
    # -------------------------
    non_gap = (aln != "-") & (aln != "N")
    column_occupancy = non_gap.sum(axis=0) / n_seq
    good_columns = column_occupancy >= min_site_occupancy

    n_good_columns = int(good_columns.sum())
    logging.info(
        f"Keeping {n_good_columns}/{n_sites} alignment columns "
        f"(min_site_occupancy={min_site_occupancy})"
    )

    if n_good_columns == 0:
        raise RuntimeError(
            "No alignment columns passed the min_site_occupancy threshold. "
            "Try lowering --min_site_occupancy."
        )

    # -------------------------
    # Step 2: compute per-sequence shared coverage
    # -------------------------
    kept = []
    removed = []
    coverage_map = {}

    denom = n_good_columns

    for i, rec in enumerate(records):
        seq_non_gap = non_gap[i]
        shared_positions = seq_non_gap & good_columns
        shared_fraction = shared_positions.sum() / denom if denom > 0 else 0.0

        coverage_map[rec.id] = shared_fraction

        if force_keep_ids and rec.id in force_keep_ids:
            kept.append(rec)
            continue

        if shared_fraction >= min_shared_fraction:
            kept.append(rec)
        else:
            removed.append(rec.id)

    logging.info(f"Kept {len(kept)} / {len(records)} sequences")

    if removed:
        logging.info("Removed sequences (low shared coverage):")
        for r in removed:
            logging.info(f"  - {r} (shared={coverage_map[r]:.3f})")

    if len(kept) == 0:
        raise RuntimeError(
            "All sequences were removed by shared coverage filtering. "
            "Try lowering --min_shared_fraction."
        )

    # -------------------------
    # Step 3: actually trim alignment columns
    # -------------------------
    trimmed_records = []

    for rec in kept:
        seq = np.array(list(str(rec.seq).upper()))
        trimmed_seq = "".join(seq[good_columns])

        rec.seq = Seq(trimmed_seq)
        trimmed_records.append(rec)

    trimmed_len = len(trimmed_records[0].seq) if trimmed_records else 0
    logging.info(f"Trimmed alignment length: {trimmed_len} bp")

    # -------------------------
    # Step 4: write filtered alignment
    # -------------------------
    with open(output_file, "w") as out:
        SeqIO.write(trimmed_records, out, "fasta")

    logging.info(f"Filtered alignment saved to: {output_file}")

    # -------------------------
    # Step 5: write reports
    # -------------------------
    removed_file = output_file.replace(".fasta", "_removed_sequences.tsv")
    with open(removed_file, "w") as f:
        f.write("sequence_id\tshared_fraction\n")
        for r in removed:
            f.write(f"{r}\t{coverage_map[r]:.6f}\n")

    logging.info(f"Removed sequences written to: {removed_file}")

    kept_file = output_file.replace(".fasta", "_kept_sequences.tsv")
    with open(kept_file, "w") as f:
        f.write("sequence_id\tshared_fraction\n")
        for r in trimmed_records:
            f.write(f"{r.id}\t{coverage_map[r.id]:.6f}\n")

    # optional: write kept columns report
    columns_file = output_file.replace(".fasta", "_kept_columns.tsv")
    with open(columns_file, "w") as f:
        f.write("original_column_index\toccupancy_fraction\tkept\n")
        for idx, occ in enumerate(column_occupancy, start=1):
            f.write(f"{idx}\t{occ:.6f}\t{int(good_columns[idx-1])}\n")

    logging.info(f"Column occupancy report written to: {columns_file}")

    return [r.id for r in trimmed_records], removed, coverage_map


# -------------------------------------------------------------------
# Calculate average neighbor identity (alignment-based, chloroplast)
# -------------------------------------------------------------------

def run_alnPairDist(alignment_file, output_file):
    alignment_file = os.path.abspath(alignment_file)
    outdir = os.path.dirname(output_file)
    os.makedirs(outdir, exist_ok=True)

    aln = AlignIO.read(alignment_file, "fasta")

    # convert once

    #records = [(rec.id, str(rec.seq).upper()) for rec in aln]
    records = [(rec.description, str(rec.seq).upper()) for rec in aln] #full description

    def calc_pid(seq1, seq2):
        matches = 0
        compared = 0

        for a, b in zip(seq1, seq2):
            # skip sites with gaps or ambiguous characters
            if a not in "ACGT" or b not in "ACGT":
                continue

            compared += 1
            if a == b:
                matches += 1

        if compared == 0:
            return 0.0

        return (matches / compared) * 100.0

    with open(output_file, "w") as out:
        # keep at least 11 columns so parse_alnPairDist() still works
        out.write("seq1\tseq2\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\tpid\n")

        for (id1, seq1), (id2, seq2) in itertools.combinations(records, 2):
            pid = calc_pid(seq1, seq2)

            out.write(
                f"{id1}\t{id2}\t.\t.\t.\t.\t.\t.\t.\t.\t{pid:.6f}\n"
            )

    logging.info(f"Saved pairwise distances to: {output_file}")
    

def assign_alignment_ani(species_dict, ani_results):
    for sp, entries in species_dict.items():
        for entry in entries:
            acc = entry.get("accession")

            if acc in ani_results:
                entry["ani"] = ani_results[acc]["ani"]

    return species_dict




def parse_alnPairDist(tsv_file, target_id, acc2species):
    results = {}

    target_species = acc2species.get(target_id, "unknown")

    with open(tsv_file) as f:
        header = next(f)

        for line in f:
            parts = line.strip().split("\t")

            if len(parts) < 11:
                continue

            s1 = parts[0]
            s2 = parts[1]
            pid = float(parts[10])
            ani = pid / 100

            # Case 1: target is s1
            if s1 == target_id:
                query = s2

            # Case 2: target is s2
            elif s2 == target_id:
                query = s1

            else:
                continue  # skip comparisons not involving target

            results[query] = {
                "ani": ani,
                "target_species": target_species,
                "query_species": acc2species.get(query, "unknown")
            }

    return results
 


def write_pairwise_with_species(dist_file, output_file, acc2species, target_id):
    target_species = acc2species.get(target_id, "unknown")

    with open(dist_file) as infile, open(output_file, "w") as out:
        header = next(infile)

        out.write("target_species\tquery_species\tANI\n")

        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue

            s1 = parts[0]
            s2 = parts[1]
            pid = float(parts[10])
            ani = pid / 100

            if s1 == target_id:
                query = s2
            elif s2 == target_id:
                query = s1
            else:
                continue

            out.write(
                f"{target_species}\t{acc2species.get(query,'unknown')}\t{ani}\n"
            )



# --------------------------------------------------
# Calculate mash distances (mitochondrial, nuclear)
# --------------------------------------------------

def compute_mash_distances(genome_files, target_acc):
    """
    genome_files:
        {acc: {"path": ..., "species": ...}}
    """
    if target_acc not in genome_files:
        raise ValueError(f"Target accession {target_acc} not found in genome_files")

    results = {}
    target_path = genome_files[target_acc]["path"]
    target_msh = mash_sketch(target_path)

    for acc, meta in genome_files.items():
        if acc == target_acc:
            continue

        msh = mash_sketch(meta["path"])
        cmd = f'mash dist "{target_msh}" "{msh}"'
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if res.returncode != 0:
            logging.warning(f"Mash failed for {target_acc} vs {acc}: {res.stderr}")
            continue

        if not res.stdout.strip():
            logging.warning(f"No Mash output for {target_acc} vs {acc}")
            continue

        fields = res.stdout.strip().split()
        if len(fields) < 3:
            logging.warning(f"Unexpected Mash output for {target_acc} vs {acc}: {res.stdout}")
            continue

        dist = float(fields[2])

        results[acc] = {
            "mash_distance": dist
        }

    return results


def compute_mash_distances_from_fasta_files(fasta_files, target_acc):
    """
    Compute Mash distances between one target FASTA and all other FASTA files 
    for mitochondrial genome sequence(s).

    fasta_files:
        {accession: fasta_path}

    Returns:
        {accession: {"mash_distance": dist}}
    """
    if target_acc not in fasta_files:
        raise ValueError(f"Target accession {target_acc} not found in fasta_files")

    results = {}
    target_path = fasta_files[target_acc]
    target_msh = mash_sketch(target_path)

    for acc, fasta_path in fasta_files.items():
        if acc == target_acc:
            continue

        msh = mash_sketch(fasta_path)
        cmd = f'mash dist "{target_msh}" "{msh}"'
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if not res.stdout.strip():
            logging.warning(f"No Mash output for {target_acc} vs {acc}")
            continue

        fields = res.stdout.strip().split()
        dist = float(fields[2])

        results[acc] = {
            "mash_distance": dist
        }

    return results



def mash_sketch(genome_path):
    """
    Create a Mash sketch for a nuclear genome.
    """
    if genome_path in msh_cache:
        return msh_cache[genome_path]

    msh_file = genome_path + ".msh"

    if not os.path.exists(msh_file):
        subprocess.run(f'mash sketch "{genome_path}"', shell=True, check=True)

    msh_cache[genome_path] = msh_file
    return msh_file



def mash_select_best(genome_files):
    """
    Select best representative genome (medoid) based on pairwise Mash distances.

    Returns:
        (best_accession, best_score)
    """
    accs = list(genome_files.keys())

    if not accs:
        raise RuntimeError("No genomes available for Mash selection.")

    if len(accs) == 1:
        return accs[0], 0.0

    msh_map = {}
    for acc, meta in genome_files.items():
        msh_map[acc] = mash_sketch(meta["path"])

    distances = {acc: [] for acc in accs}

    for i in range(len(accs)):
        for j in range(i + 1, len(accs)):
            acc1, acc2 = accs[i], accs[j]
            msh1, msh2 = msh_map[acc1], msh_map[acc2]

            cmd = f'mash dist "{msh1}" "{msh2}"'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if not result.stdout.strip():
                continue

            dist = float(result.stdout.split()[2])

            logging.info(f"Mash distance {acc1} vs {acc2}: {dist}")

            distances[acc1].append(dist)
            distances[acc2].append(dist)

    best_acc = None
    best_score = float("inf")

    for acc, dists in distances.items():
        if not dists:
            continue

        mean_dist = sum(dists) / len(dists)
        logging.info(f"{acc}: mean Mash distance = {mean_dist}")

        if mean_dist < best_score:
            best_score = mean_dist
            best_acc = acc

    if best_acc is None:
        best_acc = accs[0]
        best_score = 0.0

    logging.info(f"Selected representative genome: {best_acc}")
    return best_acc, best_score


def assign_mash_results(species_dict, mash_results):
    for sp, entries in species_dict.items():
        for entry in entries:
            acc = entry.get("accession")
            if acc in mash_results:
                entry["mash_distance"] = mash_results[acc]["mash_distance"]

    return species_dict



# -------------------------------------
# Collapse entries per species by mean
# -------------------------------------

def collapse_to_species(outfolder, taxon_name, genome_type=None):
    """
    Collapse TSV results to species level.

    Rules:
    - chloroplast: mean ANI per target/query species pair
    - mitochondrial: mean mash_distance per target/query species pair
    - nuclear_genome: mean mash_distance per target/query species pair
    """
    taxon_clean = taxon_name.replace(" ", "_")
    tsv_files = glob.glob(os.path.join(outfolder, f"{taxon_clean}_*.tsv"))
    tsv_files = [f for f in tsv_files if not f.endswith("_species_collapsed.tsv")]

    logging.info(f"collapse_to_species() called for taxon={taxon_name}, genome_type={genome_type}")
    logging.info(f"Found {len(tsv_files)} TSV file(s) for collapsing: {tsv_files}")

    if not tsv_files:
        logging.warning("No TSV files found for collapsing.")
        return

    for tsv_file in tsv_files:
        try:
            df = pd.read_csv(tsv_file, sep="\t")
        except Exception as e:
            logging.warning(f"Could not read {tsv_file}: {e}")
            continue

        col_rename_map = {
            "Target_Species": "target_species",
            "Other_Species": "query_species",
            "ANI": "ani",
            "Mash_Distance": "mash_distance",
            "Target_Seq": "target_accession",
            "Other_Seq": "query_accession",
            "Species": "query_species"
        }
        df = df.rename(columns={k: v for k, v in col_rename_map.items() if k in df.columns})

        if genome_type == "chloroplast":
            metric = "ani"
            ascending = False
        elif genome_type in ["mitochondrial", "nuclear_genome"]:
            metric = "mash_distance"
            ascending = True
        else:
            logging.warning(f"Unknown genome_type '{genome_type}', skipping {tsv_file}")
            continue

        if "target_species" not in df.columns:
            if "Target_Taxon" in df.columns:
                df["target_species"] = df["Target_Taxon"]
            else:
                df["target_species"] = taxon_name

        required_cols = {"target_species", "query_species", metric}
        missing = required_cols - set(df.columns)

        if missing:
            logging.warning(f"Skipping {tsv_file}: missing required columns {missing}")
            continue

        df = df.dropna(subset=[metric])

        if df.empty:
            logging.warning(f"Skipping {tsv_file}: no valid values in metric column '{metric}'")
            continue

        collapsed = (
            df.groupby(["target_species", "query_species"], as_index=False)[metric]
            .mean()
            .sort_values(by=metric, ascending=ascending)
        )

        #best_only = collapsed.groupby("target_species", as_index=False).first()

        out_name = tsv_file.replace(".tsv", "_species_collapsed.tsv")
        collapsed.to_csv(out_name, sep="\t", index=False)

        logging.info(f"Saved species-collapsed TSV: {out_name}")




# --------------
# Output writer
# --------------
def write_output(species_dict, out, target_id=None, taxon_name=None, genome_type=None):
    """
    Write species_dict to a TSV file.
    """
    with open(out, "w") as f:
        if genome_type == "chloroplast":
            headers = [
                "Species", "Accession", "Assembly_Name",
                "Assembly_Level", "Assembly_Date", "Provider",
                "Bioproject", "Taxonomy",
                "Sequence_Length", "Accession_Taxid", "ANI",
                "Target_Accession", "Target_Taxon"
            ]
        elif genome_type == "mitochondrial":
            headers = [
                "Species", "Accession", "Assembly_Name",
                "Assembly_Level", "Assembly_Date", "Provider",
                "Bioproject", "Taxonomy",
                "Sequence_Length", "Accession_Taxid", "Mash_Distance",
                "Target_Accession", "Target_Taxon"
            ]
        else:
            headers = [
                "Species", "Accession", "Assembly_Name", "Assembly_Level",
                "Assembly_Method", "Assembly_Date", "Provider",
                "Bioproject", "Biosample", "Annotated",
                "N50", "Genome_Size", "Chromosomes", "Mash_Distance",
                "Target_Accession", "Target_Taxon"
            ]

        f.write("\t".join(headers) + "\n")

        for sp, entries in species_dict.items():
            for entry in entries:
                if isinstance(entry, StringElement):
                    entry = str(entry)

                if genome_type == "chloroplast" and isinstance(entry, dict):
                    row = [
                        entry.get("species", ""),
                        entry.get("accession", ""),
                        entry.get("assembly_name", ""),
                        entry.get("assembly_level", ""),
                        entry.get("assembly_date", ""),
                        entry.get("provider", ""),
                        entry.get("bioproject", ""),
                        entry.get("taxonomy", ""),
                        entry.get("sequence_length", ""),
                        entry.get("accession_taxid", ""),
                        entry.get("ani", ""),
                        target_id or "",
                        taxon_name or ""
                    ]

                elif genome_type == "mitochondrial" and isinstance(entry, dict):
                    row = [
                        entry.get("species", ""),
                        entry.get("accession", ""),
                        entry.get("assembly_name", ""),
                        entry.get("assembly_level", ""),
                        entry.get("assembly_date", ""),
                        entry.get("provider", ""),
                        entry.get("bioproject", ""),
                        entry.get("taxonomy", ""),
                        entry.get("sequence_length", ""),
                        entry.get("accession_taxid", ""),
                        entry.get("mash_distance", ""),
                        target_id or "",
                        taxon_name or ""
                    ]

                elif isinstance(entry, dict):
                    row = [
                        entry.get("species", ""),
                        entry.get("accession", ""),
                        entry.get("assembly_name", ""),
                        entry.get("assembly_level", ""),
                        entry.get("assembly_method", ""),
                        entry.get("assembly_date", ""),
                        entry.get("provider", ""),
                        entry.get("bioproject", ""),
                        entry.get("biosample", ""),
                        entry.get("annotation_available", ""),
                        entry.get("n50", ""),
                        entry.get("genome_size", ""),
                        entry.get("chromosomes", ""),
                        entry.get("mash_distance", ""),
                        target_id or "",
                        taxon_name or ""
                    ]

                elif isinstance(entry, str):
                    row = [sp, entry] + [""] * (len(headers) - 2)

                else:
                    continue

                f.write("\t".join(str(x) if x is not None else "" for x in row) + "\n")


def safe_int(value, default=0):
    try:
        if value is None:
            return default
        return int(str(value).replace(",", "").strip())
    except Exception:
        return default


def get_genus_from_species_name(name):
    parts = str(name).strip().split()
    return parts[0] if parts else ""


def flatten_species_dict(species_dict):
    genomes = []
    for entries in species_dict.values():
        genomes.extend(entries)
    return genomes


def select_largest_accession_from_metadata(species_dict, genome_type, preferred_genus=None):
    """
    Select the accession with the largest genome/assembly span.

    For organellar genomes: uses sequence_length.
    For nuclear genomes: uses genome_size.
    If preferred_genus is given, restrict to that genus first.
    """
    all_entries = flatten_species_dict(species_dict)

    if preferred_genus:
        preferred_genus = preferred_genus.strip().lower()
        genus_entries = [
            e for e in all_entries
            if get_genus_from_species_name(e.get("species", "")).lower() == preferred_genus
        ]
        if genus_entries:
            all_entries = genus_entries

    if genome_type in ["chloroplast", "mitochondrial"]:
        size_key = "sequence_length"
    else:
        size_key = "genome_size"

    best_entry = None
    best_size = -1

    for e in all_entries:
        size = safe_int(e.get(size_key, 0), default=0)
        if size > best_size and e.get("accession"):
            best_size = size
            best_entry = e

    if best_entry is None:
        raise RuntimeError("Could not select target accession from metadata.")

    logging.info(
        f"Selected largest target accession from metadata: {best_entry['accession']} "
        f"(species={best_entry.get('species','unknown')}, size={best_size})"
    )
    return best_entry["accession"]

def restrict_species_dict_to_accessions(species_dict, accessions):
    accessions = set(accessions)
    filtered = {}

    for species, entries in species_dict.items():
        kept = [e for e in entries if e.get("accession") in accessions]
        if kept:
            filtered[species] = kept

    return filtered


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def main():
    parser = ArgumentParser(description="Find the closest available reference genome sequence(s) of a given taxon in NCBI.")
    parser.add_argument('--taxon', required=True, help="Species or higher-level taxon name (e.g., Genus or Family).")
    parser.add_argument('--outfolder', required=True, help="Output folder for the result file.")
    parser.add_argument('--genome_type', required=True,
                        choices=['chloroplast', 'mitochondrial', 'nuclear_genome'])
    parser.add_argument('--annotated', action='store_true',
                        help="Select only gene-annotated nuclear genomes.")
    parser.add_argument('--assembly_level', required=False,
                        choices=['contig', 'scaffold', 'chromosome', 'complete'],
                        help="Choose the assembly level of the nuclear genome (STRING; contig scaffold, chromosome, complete).")
    parser.add_argument("--genome_size_mb", type=float, default=None,
                        help=("Expected genome size in Mb. "
                              "If the target taxon has no genome, this value is used "
                              "to select the most similar proxy genome (FLOAT, default: None)."))
    parser.add_argument("--max_genomes", required=False, type=int, default=None,
                        help="Maximum number of genomes to compare target genome."
                        "The more genomes allowed, the longer the distance-based sample calculations will take"
                        "(INTEGER, default: None).")
    parser.add_argument("--min_shared_fraction", type=float, default=0.70, 
                        help="Only plastid genome sequences: The fraction of positions where sequence overlaps "
                        "with others (FLOAT, default: 0.7).")
    parser.add_argument("--min_site_occupancy", type=float, default=0.50, 
                        help="Only plastid genome sequences: The fraction of sequences that "
                        "must have a base at a site (FLOAT, default: 0.5).")
    parser.add_argument("--collapse_to_species", action="store_true", 
                        help="Collapse multiple accessions per species "
                        "to the best Average Nucleotide Identity (ANI) or Mash Distance hit.")
    parser.add_argument("--include", required=False, default='genome',
                        choices=['genome', 'cds'],
                        help="Download nuclear genome sequence(s) from NCBI (STRING; default=genome).")
    parser.add_argument("--rank", default="genus", choices=["species", "genus", "family", "order"],
                        help=("Highest taxonomic rank to include when searching for related genomes. "
                              "For species input, search expands progressively up to this rank "
                              "(STRING; default: genus)."))
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--email', required=True)
    args = parser.parse_args()
    

    validate_inputs(args.taxon, args.genome_type)
    setup_folder(args.outfolder, args.overwrite)

    tmp_dir = get_tmpdir(args.outfolder)

    log_file = os.path.join(args.outfolder, "run.log")
    logging.getLogger().handlers.clear()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )

    try:
        rank, taxid = get_taxid_and_rank(args.taxon, args.email)
    except Exception as e:
        logging.error(f"Could not resolve taxon '{args.taxon}': {e}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return

    logging.info(f"Taxon '{args.taxon}' resolved to rank='{rank}', taxid={taxid}")

    taxon_clean = args.taxon.replace(" ", "_")
    out = os.path.join(args.outfolder, f"{taxon_clean}_{args.genome_type}.tsv")

    species_dict = {}
    target_genome = None
    target_id = None

    # ---------------------------------------------------------
    # 1. Collect candidate genomes exactly once
    # ---------------------------------------------------------
    if rank == "species":
        logging.info("Species input detected")

        allowed_ranks = get_search_ranks_upto(args.rank)
        search_ranks = trim_ranks_from_input_rank("species", allowed_ranks)

        logging.info(f"User-selected maximum rank: {args.rank}")
        logging.info(f"Effective search ranks: {search_ranks}")

        species_dict = collect_genomes_across_ranks(
            taxon_name=args.taxon,
            taxid=taxid,
            genome_type=args.genome_type,
            email=args.email,
            annotated=args.annotated,
            assembly_level=args.assembly_level,
            ranks=search_ranks,
            max_genomes=args.max_genomes
        )

        if not species_dict:
            logging.error(f"No genomes found for taxon '{args.taxon}' up to rank '{args.rank}'")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return

        all_genomes = []
        for genomes in species_dict.values():
            all_genomes.extend(genomes)

        target_name = args.taxon.strip().lower()
        exact_matches = [
            g for g in all_genomes
            if g.get("species", "").strip().lower() == target_name
        ]

        if exact_matches:
            target_genome = select_largest_entry(exact_matches, args.genome_type)
            target_id = target_genome["accession"]
            logging.info(f"Using exact-species target genome: {target_id}")
        else:
            genus_name = args.taxon.split()[0]
            genus_matches = [
                g for g in all_genomes
                if g.get("species", "").lower().startswith(genus_name.lower() + " ")
            ]

            if genus_matches:
                target_genome = select_largest_entry(genus_matches, args.genome_type)
                target_id = target_genome["accession"]
                logging.info(f"No exact species genome found -> using largest genome from genus '{genus_name}': {target_id}")
            else:
                target_genome = select_largest_entry(all_genomes, args.genome_type)
                target_id = target_genome["accession"]
                logging.info(f"No genus genome found -> using largest available genome overall: {target_id}")

    else:
        logging.info("Higher-level taxon detected")

        if args.genome_type in ["chloroplast", "mitochondrial"]:
            species_dict = list_entire_taxon_genomes(taxid, args.genome_type, args.email)
        else:
            species_dict = list_taxon_nuclear_genomes(
                taxid,
                rank,
                annotated=args.annotated,
                assembly_level=args.assembly_level
            )

        if not species_dict:
            logging.error(f"No genomes found for higher-level taxon '{args.taxon}'")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return

        all_genomes = []
        for genomes in species_dict.values():
            all_genomes.extend(genomes)

        target_genome = select_largest_entry(all_genomes, args.genome_type)
        target_id = target_genome["accession"]
        logging.info(f"Using largest available genome within '{args.taxon}' as target: {target_id}")

    # ---------------------------------------------------------
    # 2. Restrict to target + remaining genomes up to max_genomes
    # ---------------------------------------------------------
    all_genomes = []
    for genomes in species_dict.values():
        all_genomes.extend(genomes)

    remaining_genomes = [g for g in all_genomes if g["accession"] != target_id]

    if args.max_genomes is None:
        selected_genomes = [target_genome] + remaining_genomes
    else:
        selected_genomes = [target_genome] + remaining_genomes[:max(0, args.max_genomes - 1)]

    accessions = [g["accession"] for g in selected_genomes]
    logging.info(f"Processing {len(accessions)} genomes")

    # Optional: reduce species_dict to selected accessions only
    filtered_species_dict = {}
    selected_set = set(accessions)
    for sp, entries in species_dict.items():
        kept = [e for e in entries if e.get("accession") in selected_set]
        if kept:
            filtered_species_dict[sp] = kept
    species_dict = filtered_species_dict

    # ---------------------------------------------------------
    # 3. Genome-type-specific analysis
    # ---------------------------------------------------------
    if args.genome_type == "nuclear_genome":
        genome_files = download_nuclear_genomes_batch(
            accessions,
            species_dict,
            tmp_dir,
            dehydrated=True,
            chunk_size=3,
            retries=4,
            include=args.include
        )

        if not genome_files:
            logging.error("No nuclear genomes downloaded successfully.")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return

        if target_id not in genome_files:
            logging.warning(f"Target accession {target_id} not downloaded -> selecting fallback")
            target_id = select_target_genome(
                genome_files,
                args.taxon,
                genome_size_mb=args.genome_size_mb
            )

        mash_results = compute_mash_distances(genome_files, target_id)
        species_dict = assign_mash_results(species_dict, mash_results)

        # flatten entries
        all_entries = []
        for sp, entries in species_dict.items():
            for e in entries:
                all_entries.append(e)

        # sort by mash distance (ascending)
        all_entries_sorted = sorted(
            all_entries,
            key=lambda e: e.get("mash_distance", float("inf"))
        )

        sorted_species_dict = {}
        for e in all_entries_sorted:
            sp = e["species"]
            sorted_species_dict.setdefault(sp, []).append(e)

        acc2species = build_accession_to_species_map(species_dict)
        actual_target_taxon = acc2species.get(target_id, "unknown")

        write_output(
            sorted_species_dict,
            out,
            target_id=target_id,
            taxon_name=actual_target_taxon,
            genome_type=args.genome_type
        )


    elif args.genome_type == "mitochondrial":
        fasta_files = download_fasta_parallel(
            accessions,
            tmp_dir,
            args.email,
            max_workers=6)

        if not fasta_files:
            logging.error("No mitochondrial FASTA files downloaded successfully.")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return

        if target_id not in fasta_files:
            logging.warning(f"Target genome {target_id} missing -> fallback selection")
            target_id = determine_target_id(fasta_files, genome_size_mb=args.genome_size_mb)

        mash_results = compute_mash_distances_from_fasta_files(fasta_files, target_id)
        species_dict = assign_mash_results(species_dict, mash_results)

        sorted_species = dict(sorted(
            species_dict.items(),
            key=lambda x: min(e.get("mash_distance", float("inf")) for e in x[1])
        ))



        acc2species = build_accession_to_species_map(species_dict)
        actual_target_taxon = acc2species.get(target_id, "unknown")

        write_output(
            sorted_species,
            out,
            target_id=target_id,
            taxon_name=actual_target_taxon,
            genome_type=args.genome_type
        )


    else:
        fasta_files = download_fasta_parallel(
            accessions,
            tmp_dir,
            args.email,
            max_workers=6)

        if not fasta_files:
            logging.error("No chloroplast FASTA files downloaded successfully.")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return

        if target_id not in fasta_files:
            logging.warning(f"Target genome {target_id} missing -> fallback selection")
            target_id = determine_target_id(fasta_files, genome_size_mb=args.genome_size_mb)

        alignment_file = os.path.join(tmp_dir, "alignment.fasta")
        run_mafft(fasta_files, alignment_file, taxon_name=args.taxon, outfolder=args.outfolder)

        filtered_alignment = os.path.join(args.outfolder, f"{taxon_clean}_filtered_alignment.fasta")
        kept_ids, removed_ids, coverage_map = filter_alignment_by_shared_coverage(
            alignment_file,
            filtered_alignment,
            min_shared_fraction=args.min_shared_fraction,
            min_site_occupancy=args.min_site_occupancy,
            force_keep_ids={target_id}
        )

        if len(kept_ids) < 2:
            raise RuntimeError("Too few sequences remain after filtering. Try lowering thresholds.")

        dist_file = os.path.join(args.outfolder, f"{taxon_clean}_pairwise_distances.tsv")
        run_alnPairDist(filtered_alignment, dist_file)

        acc2species = build_accession_to_species_map(species_dict)
        ani_results = parse_alnPairDist(dist_file, target_id, acc2species)

        sorted_results = sorted(ani_results.items(), key=lambda x: x[1]["ani"], reverse=True)
        sorted_file = os.path.join(args.outfolder, f"{taxon_clean}_sorted_by_ani.tsv")

        with open(sorted_file, "w") as f:
            f.write("Target_Seq\tTarget_Species\tOther_Seq\tOther_Species\tANI\n")
            for acc, data in sorted_results:
                f.write(f"{target_id}\t{data.get('target_species','unknown')}\t{acc}\t{data.get('query_species','unknown')}\t{data['ani']}\n")

        pretty_dist_file = os.path.join(args.outfolder, f"{taxon_clean}_pairwise_with_species.tsv")
        write_pairwise_with_species(dist_file, pretty_dist_file, acc2species, target_id)

        species_dict = assign_alignment_ani(species_dict, ani_results)

        sorted_species = dict(sorted(
            species_dict.items(),
            key=lambda x: max(e.get("ani", 0) for e in x[1]),
            reverse=True
        ))

        write_output(
            sorted_species,
            out,
            target_id=target_id,
            taxon_name=args.taxon,
            genome_type=args.genome_type
        )

    # ---------------------------------------------------------
    # 4. Optional collapse
    # ---------------------------------------------------------
    if args.collapse_to_species:
        logging.info(f"Collapsing results to species level for genome_type={args.genome_type}")
        collapse_to_species(
            outfolder=args.outfolder,
            taxon_name=args.taxon,
            genome_type=args.genome_type
        )

    shutil.rmtree(tmp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()
