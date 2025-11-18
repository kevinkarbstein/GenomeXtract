#!/usr/bin/env python

"""
This script assembles organellar genomes from a given group (e.g., genus, family, or order) derived from the NCBI database. 

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

import argparse
import os
import subprocess
import Bio 
from Bio import SeqIO, SeqFeature, SeqRecord
import pandas as pd
import urllib.request
import zipfile
import requests
import sys 
import io 
import re


def process_genbank_files(file_paths, select_group=None):
    """
    Process the GenBank files to extract gene information.
    """
    organism_section_gene_map = {}
    all_genes = set()

    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        with open(file_path, 'r') as file:
            organism_name = "Unknown Organism"
            section_name = "Unknown Section"
            genes = {}
            current_gene = None
            in_cds_feature = False
            sequence_lines = []
            
            for line in file:
                if line.startswith("  ORGANISM  "):
                    try:
                        organism_name = " ".join(line.split()[1:3])  # Expect two words for organism name
                    except IndexError:
                        organism_name = "Unknown Organism"

                elif line.startswith("            "):
                    terms = line.split(";")
                    section_name = terms[-2].strip() if len(terms) > 1 else "Unknown Section"
                    if select_group and select_group != section_name:
                        continue  # Skip this file if it doesn't match the selected group

                if line.startswith("     CDS"):
                    in_cds_feature = True
                elif in_cds_feature and "/gene=" in line:
                    current_gene = line.split("=")[1].strip().replace('"', '').lower()
                    genes[current_gene] = None
                    all_genes.add(current_gene)
                elif in_cds_feature and line.startswith("ORIGIN"):
                    in_cds_feature = False
                    sequence_lines = []  # Start capturing sequence lines
                elif sequence_lines is not None:
                    if line.startswith("//"):
                        if current_gene and sequence_lines:
                            genes[current_gene] = ''.join(sequence_lines).replace(' ', '').replace('\n', '')
                        sequence_lines = None  # Stop sequence capture
                    else:
                        sequence_lines.append(''.join(line.strip().split()[1:]))

            organism_section_gene_map[organism_name] = {
                "section": section_name,
                "genes": genes
            }

    return organism_section_gene_map, all_genes


def reorder_columns(df, section_order):
    """
    Reorder the columns of a DataFrame based on a specified section order.
    """
    ordered_columns = ['Gene']
    for section in section_order:
        ordered_columns.extend([col for col in df.columns if col != 'Gene' and section in col])
    return df[ordered_columns]


def write_summary_files(organism_section_gene_map, all_genes, output_base, feature_section_summary, section_order, overwrite):
    """
    Write the summary information to CSV and XLSX files.
    """
    data = {'Gene': sorted(all_genes)}
    for organism in organism_section_gene_map:
        data[organism] = [gene if gene in organism_section_gene_map[organism]['genes'] else '' for gene in sorted(all_genes)]
    
    df = pd.DataFrame(data)
    
    if section_order:
        df = reorder_columns(df, section_order)
    
    # Write to CSV
    csv_output_file = f"{output_base}.csv"
    if not os.path.exists(csv_output_file) or overwrite:
        df.to_csv(csv_output_file, index=False)
        print(f"Gene summary written to {csv_output_file}")
    
    # Write to XLSX
    xlsx_output_file = f"{output_base}.xlsx"
    if not os.path.exists(xlsx_output_file) or overwrite:
        df.to_excel(xlsx_output_file, index=False)
        print(f"Gene summary written to {xlsx_output_file}")
    
    if feature_section_summary:
        section_gene_map = {}
        for organism in organism_section_gene_map:
            section = organism_section_gene_map[organism]['section']
            if section not in section_gene_map:
                section_gene_map[section] = set()
            section_gene_map[section].update(organism_section_gene_map[organism]['genes'].keys())
        
        section_data = {'Gene': sorted(all_genes)}
        for section in section_gene_map:
            section_data[section] = [gene if gene in section_gene_map[section] else '' for gene in sorted(all_genes)]
        
        section_df = pd.DataFrame(section_data)
        
        if section_order:
            section_df = reorder_columns(section_df, section_order)
        
        # Write section summary to CSV
        section_csv_output_file = f"{output_base}_section_summary.csv"
        if not os.path.exists(section_csv_output_file) or overwrite:
            section_df.to_csv(section_csv_output_file, index=False)
            print(f"Section-wise gene summary written to {section_csv_output_file}")
        
        # Write section summary to XLSX
        section_xlsx_output_file = f"{output_base}_section_summary.xlsx"
        if not os.path.exists(section_xlsx_output_file) or overwrite:
            section_df.to_excel(section_xlsx_output_file, index=False)
            print(f"Section-wise gene summary written to {section_xlsx_output_file}")


def write_combined_gene_sequences(organism_section_gene_map, output_fasta_file, overwrite):
    """
    Write all gene sequences from all organisms into a single combined FASTA file.
    """
    output_dir = os.path.dirname(output_fasta_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if not os.path.exists(output_fasta_file) or overwrite:
        with open(output_fasta_file, 'w') as f:
            for organism, info in organism_section_gene_map.items():
                for gene, sequence in info['genes'].items():
                    if sequence:
                        f.write(f">{organism.replace(' ', '_')}|{gene}\n{sequence}\n")
        print(f"Combined gene sequences written to {output_fasta_file}")
    else:
        print(f"Combined FASTA file already exists: {output_fasta_file}. Use --overwrite to regenerate.")


def align_combined_sequences(input_fasta_file, aligned_output_file, overwrite):
    """
    Align combined gene sequences using MAFFT.
    """
    if not os.path.exists(aligned_output_file) or overwrite:
        mafft_cmd = f"mafft --thread -1 --auto {input_fasta_file} > {aligned_output_file}"
        subprocess.run(mafft_cmd, shell=True)
        print(f"Aligned sequences written to {aligned_output_file}")


def run_raxml_ng(aligned_fasta_file, output_dir, overwrite, raxml_model="GTR+G", bootstrap_replicates=100):
    """
    Run RAxML-NG for phylogenetic analysis on the aligned sequences.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, "raxml_output")
    if not os.path.exists(output_path) or overwrite:
        raxml_cmd = f"raxml-ng --all --msa {aligned_fasta_file} --model {raxml_model} --bs-trees {bootstrap_replicates} --prefix {output_path} --threads AUTO --bs-metric fbp,tbe"
        subprocess.run(raxml_cmd, shell=True)
        print(f"RAxML-NG analysis completed for {aligned_fasta_file}")

def run_iqtree(aligned_fasta_file, output_dir, overwrite, bootstrap_replicates=100):
    """
    Run IQ-TREE for phylogenetic analysis on the aligned sequences.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, "iqtree_output")
    if not os.path.exists(output_path) or overwrite:
        iqtree_cmd = f"iqtree -s {aligned_fasta_file} -pre concat -T AUTO -b {bootstrap_replicates} --tbe --prefix {output_path}"
        subprocess.run(iqtree_cmd, shell=True)
        print(f"IQ-TREE analysis completed for {aligned_fasta_file}")



def main():
    parser = argparse.ArgumentParser(description="Process GenBank files to extract gene information and perform phylogenetic analyses.")
    parser.add_argument('--input', nargs='+', required=True, help='Path to the GenBank files')
    parser.add_argument('--feature_summary', required=True, help='File name for feature summary CSV and XLSX files')
    parser.add_argument('--group_feature_summary', action='store_true', help='Whether to generate a section-wise feature summary')
    parser.add_argument('--group_order', nargs='+', help='Optional: Order of sections in the summary files')
    parser.add_argument('--generate_gene_sequences', action='store_true', help='Generate gene sequences in FASTA format')
    parser.add_argument('--align_sequences', action='store_true', help='Align gene sequences using MAFFT')
    parser.add_argument('--run_raxml', action='store_true', help='Run RAxML-NG for phylogenetic analysis')
    parser.add_argument('--raxml_model', default='GTR+G', help='RAxML-NG model to use (default: GTR+G)')
    parser.add_argument('--bootstrap_replicates', type=int, default=100, help='Number of bootstrap replicates for RAxML-NG (default: 100)')
    parser.add_argument('--run_iqtree', action='store_true', help='Run IQ-TREE for phylogenetic analysis')
    parser.add_argument('--select_group', default=None, help='Limit gene extraction to a specific group')
    parser.add_argument('--output_dir', default='gene_output', help='Output directory for gene sequences and alignment results')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files if they exist')
    
    args = parser.parse_args()

    # Process GenBank files and extract gene information
    organism_section_gene_map, all_genes = process_genbank_files(args.input, select_group=args.select_group)
    
    # Write summary files
    write_summary_files(organism_section_gene_map, all_genes, args.feature_summary, args.group_feature_summary, args.group_order, args.overwrite)
    
    # Generate combined gene sequence FASTA file
    combined_fasta_file = os.path.join(args.output_dir, "combined_gene_sequences.fasta")
    write_combined_gene_sequences(organism_section_gene_map, combined_fasta_file, args.overwrite)

    # Align combined sequences if specified
    aligned_fasta_file = os.path.join(args.output_dir, "aligned_gene_sequences.fasta")
    if args.align_sequences:
        align_combined_sequences(combined_fasta_file, aligned_fasta_file, args.overwrite)
    
    # Run RAxML-NG if specified
    if args.run_raxml:
        run_raxml_ng(aligned_fasta_file, args.output_dir, args.overwrite, args.raxml_model, args.bootstrap_replicates)
    # Run IQ-TREE if specified
    if args.run_iqtree:
        run_iqtree(aligned_fasta_file, args.output_dir, args.overwrite, args.bootstrap_replicates)


if __name__ == "__main__":
    main()
