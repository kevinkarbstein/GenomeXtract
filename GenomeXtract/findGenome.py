#!/usr/bin/env python

"""
This script finds all available plastomes, mitogenomes, and nuclear genomes from a given group (e.g., genus, family, or order) in the public NCBI database. 

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
import logging
import zipfile
import json
import math
import shutil
import subprocess
from argparse import ArgumentParser
from Bio import Entrez, SeqIO


# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def validate_inputs(group, genome_type):
    if not group or not genome_type:
        raise ValueError("Both group and genome_type must be provided.")

def setup_folder(outfolder, overwrite):
    if os.path.exists(outfolder):
        if overwrite:
            logging.info(f"Overwriting existing folder: {outfolder}")
            shutil.rmtree(outfolder)
        else:
            logging.warning(f"Folder {outfolder} already exists. Use --overwrite to replace it.")
            raise FileExistsError(f"Folder {outfolder} exists.")
    os.makedirs(outfolder, exist_ok=True)

def setup_excluded_folder(outfolder):
    excluded_folder = os.path.join(outfolder, "excluded")
    os.makedirs(excluded_folder, exist_ok=True)
    return excluded_folder

def estimate_organellar_genome_size(term, retmax, db):
    """Estimate the total size of the genomes in gigabytes."""
    try:
        handle = Entrez.esearch(db=db, term=term, idtype="acc", retmax=retmax)
        record = Entrez.read(handle)
        id_list = record['IdList']
        total_size = 0
        for genome_id in id_list:
            handle = Entrez.esummary(db=db, id=genome_id, rettype="gb", retmode="text")
            summary = Entrez.read(handle)[0]
            total_size += int(summary['Length'])
        total_size_gb = total_size / (1024 ** 3)  # Convert from bytes to gigabytes
        return total_size_gb
    except Exception as e:
        logging.error(f"Error estimating genome size: {e}")
        raise

def remove_duplicates(outfolder):
    """Remove .gb files with duplicate FASTA sequences, prioritizing NC_* or the latest release."""
    fasta_sequences = {}
    excluded_folder = setup_excluded_folder(outfolder)
    
    for filename in os.listdir(outfolder):
        if filename.endswith(".gb"):
            filepath = os.path.join(outfolder, filename)
            with open(filepath, "r") as file:
                record = SeqIO.read(file, "genbank")
                seq = str(record.seq)
                organism = record.annotations.get("organism", "Unknown")

                if seq in fasta_sequences:
                    fasta_sequences[seq].append((filepath, organism))
                else:
                    fasta_sequences[seq] = [(filepath, organism)]
    
    for seq, files in fasta_sequences.items():
        if len(files) > 1:
            files.sort(key=lambda f: (
                not os.path.basename(f[0]).startswith("NC_"),
                SeqIO.read(open(f[0]), "genbank").annotations.get("date", "")
            ))
            kept_file = files[0]
            for file_to_remove in files[1:]:
                shutil.move(file_to_remove[0], os.path.join(excluded_folder, os.path.basename(file_to_remove[0])))
                logging.info(f"Removed duplicate '{file_to_remove[0]}', kept '{kept_file[0]}'.")

def apply_max_individuals_per_species(outfolder, max_individuals_per_species):
    """Keep only the most recent files if max individuals per species is specified."""
    organisms = {}
    excluded_folder = setup_excluded_folder(outfolder)
    
    for filename in os.listdir(outfolder):
        if filename.endswith(".gb"):
            filepath = os.path.join(outfolder, filename)
            with open(filepath, "r") as file:
                record = SeqIO.read(file, "genbank")
                organism = record.annotations.get("organism", "Unknown")

                if organism not in organisms:
                    organisms[organism] = []
                organisms[organism].append(filepath)
    
    for organism, files in organisms.items():
        if len(files) > max_individuals_per_species:
            files.sort(key=lambda f: SeqIO.read(open(f), "genbank").annotations.get("date", ""))
            for file_to_remove in files[:-max_individuals_per_species]:
                shutil.move(file_to_remove, os.path.join(excluded_folder, os.path.basename(file_to_remove)))
                logging.info(f"Removed '{file_to_remove}' due to max individuals limit for '{organism}'.")

def fetch_genome_ids(term, retmax, batch_size, db):
    """Fetch genome IDs from NCBI in batches."""
    try:
        handle = Entrez.esearch(db=db, term=term, idtype="acc", retmax=retmax)
        record = Entrez.read(handle)
        return record['IdList'], math.ceil(len(record['IdList']) / batch_size)
    except Exception as e:
        logging.error(f"Error fetching genome IDs: {e}")
        raise

def download_genomes(batch_ids, db):
    """Download genomes from NCBI."""
    try:
        handle = Entrez.efetch(db=db, id=batch_ids, rettype="gb", retmode="text")
        return SeqIO.parse(handle, "genbank")
    except Exception as e:
        logging.error(f"Error downloading genomes: {e}")
        raise

def save_genomes(records, outfolder, genome_type):
    """Save downloaded genomes to files."""
    for record in records:
        filename = os.path.join(outfolder, f"{record.name}_{genome_type}.gb")
        if os.path.exists(filename):
            logging.warning(f"File '{filename}' already exists. Skipping.")
            continue
        with open(filename, "w") as f:
            SeqIO.write(record, f, "genbank")
        logging.info(f"Saved {record.name} to '{filename}'.")

def find_full_organelle_genome(group, outfolder, genome_type, duplicate_removal, max_individuals_per_species, overwrite, email):
    """
    Download and process organelle genomes (chloroplast or mitochondrial).
    """
    Entrez.email = email  # Set email for NCBI Entrez
    
    # Construct search term based on genome type
    db = 'nucleotide'
    term = f'("{group}"[Organism]) AND ("complete genome")'
    if genome_type == 'chloroplast':
        term += ' AND chloroplast'
    elif genome_type == 'mitochondrial':
        term += ' AND mitochondrion[Title]'
    else:
        logging.error(f"Unsupported genome type: {genome_type}")
        raise ValueError("Invalid genome_type. Use 'chloroplast', 'mitochondrial', or 'nuclear'.")

    try:
        # Estimate download size
        estimated_size_gb = estimate_organellar_genome_size(term=term, retmax=100000, db=db)
        logging.info(f"Estimated download size: {estimated_size_gb:.2f} GB")
        
        #user_input = input(f"Do you want to proceed? Estimated size: {estimated_size_gb:.2f} GB (yes/no): ").strip().lower()
        #if user_input != 'yes':
            #logging.info("Download canceled by user.")
            #return

        # Fetch genome IDs and batches
        ids, num_batches = fetch_genome_ids(term=term, retmax=100000, batch_size=50, db=db)
        logging.info(f"Found {len(ids)} genomes for group '{group}' and type '{genome_type}'.")

        # Process batches of genomes
        for i in range(num_batches):
            batch_ids = ids[i * 100:(i + 1) * 100]
            logging.info(f"Downloading batch {i + 1}/{num_batches} with {len(batch_ids)} IDs.")
            records = download_genomes(batch_ids=batch_ids, db=db)
            save_genomes(records, outfolder=outfolder, genome_type=genome_type)

        # Remove duplicates if enabled
        if duplicate_removal:
            remove_duplicates(outfolder)

        # Apply max individuals filter if specified
        if max_individuals_per_species:
            apply_max_individuals_per_species(outfolder, max_individuals_per_species)

        logging.info(f"Download and processing complete. Genomes saved to {outfolder}.")
    
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise


def preview_nuclear_genome_size(group, assembly_level):
    """Preview the nuclear genome dataset size using the NCBI Datasets API."""
    try:
        command = f"""
        datasets download genome taxon "{group}" --preview --include genome --annotated --assembly-level "{assembly_level}" --assembly-version latest --exclude-atypical
        """
        logging.info(f"Running preview command: {command.strip()}")
        result = subprocess.run(command, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            # Extract JSON output and parse it
            logging.info(f"Preview command output: {result.stdout}")
            for line in result.stdout.splitlines():
                if line.strip().startswith('{') and line.strip().endswith('}'):
                    json_output = json.loads(line.strip())
                    estimated_size_mb = json_output.get("estimated_file_size_mb", 0)
                    return estimated_size_mb / 1024  # Convert MB to GB
        else:
            logging.error(f"Error in preview command. Return code: {result.returncode}. Stderr: {result.stderr}")
            return None

    except Exception as e:
        logging.error(f"Failed to preview nuclear genome size: {str(e)}")
        return None

def find_nuclear_genomes(group, outfolder, genome_type, annotated, assembly_level):
    """Search, download, and process genomes using NCBI Datasets API."""
    if genome_type == "nuclear_genome":
        logging.info(f"Handling 'nuclear' genome type for group '{group}'.")

        # Preview genome size
        estimated_size_gb = preview_nuclear_genome_size(group=group, assembly_level=assembly_level)
        if estimated_size_gb:
            logging.info(f"Estimated size of the single (compressed) nuclear genome fasta file: {estimated_size_gb:.2f} GB")
            #user_input = input(f"Do you want to proceed? Estimated size: {estimated_size_gb:.2f} GB (yes/no): ").strip().lower()
            #if user_input != 'yes':
                #logging.info("Download canceled by user.")
                #return

        # Handle nuclear genomes
        genomes_folder = os.path.join(outfolder, f"{genome_type}s")
        os.makedirs(genomes_folder, exist_ok=True)
        
        # Select annotated nuclear genomes
        annotated_arg = "--annotated" if annotated else ""

        zip_filename = os.path.join(outfolder, f"{group}_nuclear_genome.zip")
        command = f"""
        datasets download genome taxon "{group}" --include genome,gbff {annotated_arg} --assembly-level "{assembly_level}" --assembly-version latest --exclude-atypical --filename {zip_filename}
        """

        try:
            logging.info(f"Running NCBI Datasets command: {command.strip()}")
            result = subprocess.run(command, shell=True, executable='/bin/bash')

            if result.returncode == 0:
                logging.info(f"Nuclear genomes for group '{group}' downloaded successfully.")

                # Unzip the file
                logging.info(f"Unzipping {zip_filename}...")
                with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                    zip_ref.extractall(genomes_folder)
                logging.info(f"Successfully unzipped {zip_filename}.")

                # Collect .fna file paths for comparison
                fna_files = [os.path.join(root, file) for root, dirs, files in os.walk(genomes_folder) for file in files if file.endswith(".fna")]

                # Remove duplicates if needed (NCBI Genome already applies strict filters)
                # if duplicate_removal:
                #     remove_duplicates(fna_files)

                # Apply max individuals per species if needed (currently, there even only a few genomes within a plant order)
                # if max_individuals_per_species:
                #     apply_max_individuals_per_species(fna_files, max_individuals_per_species)

                return genomes_folder
                
            else:
                logging.error(f"Error downloading genomes. Return code: {result.returncode}")

        except Exception as e:
            logging.error(f"Failed to download or unzip nuclear genomes: {str(e)}")

    return None  # Return None if the process fails or isn't for nuclear genomes


def main():
    parser = ArgumentParser(description="Script for downloading and processing genomes.")
    parser.add_argument('-g', '--group', required=True, help="Taxonomic group or organism name.")
    parser.add_argument('-o', '--outfolder', required=True, help="Output folder for downloaded files.")
    parser.add_argument('-t', '--genome_type', required=True, choices=['chloroplast', 'mitochondrial', 'nuclear_genome'], help="Type of genome to process.")
    parser.add_argument('--annotated', action='store_true', help="Select only gene-annotated nuclear genomes.")
    parser.add_argument('--assembly_level', required=False, choices=['scaffold', 'chromosome'], help="Choose the assmbly level of the nuclear genome.")
    parser.add_argument("--batch_size", type=int, default=50, help="Batch size for downloading.")
    parser.add_argument('--duplicate_removal', action='store_true', help="Remove duplicate files. Prioritize NC_* or the latest release.")
    parser.add_argument('--max_individuals', type=int, help="Maximum individuals per species.")
    parser.add_argument('--overwrite', action='store_true', help="Overwrite existing output folder.")
    parser.add_argument('--email', required=True, help="Your email for NCBI Entrez queries.")
    args = parser.parse_args()

    # Validate inputs
    validate_inputs(args.group, args.genome_type)

    # Set up the output folder
    setup_folder(args.outfolder, overwrite=args.overwrite)

    # Process genomes based on genome type
    try:
        if args.genome_type == "nuclear_genome":
            find_nuclear_genomes(
                group=args.group,
                annotated=args.annotated,
                assembly_level=args.assembly_level,
                outfolder=args.outfolder,
                genome_type=args.genome_type
            )
        else:
            find_full_organelle_genome(
                group=args.group,
                outfolder=args.outfolder,
                genome_type=args.genome_type,
                duplicate_removal=args.duplicate_removal,
                max_individuals_per_species=args.max_individuals,
                overwrite=args.overwrite,
                email=args.email
            )
    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
