#!/usr/bin/env python

"""
This script assembles organellar genes from a given group (e.g., genus, family, or order) derived from the NCBI database.

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
from Bio import SeqIO, SeqFeature, SeqRecord
import pandas as pd
import os
import subprocess
import shutil

def extract_genes_with_cds(file_path):
    """
    Extract gene names and sequences associated with 'CDS' features from a GenBank file.
    
    Args:
    - file_path (str): Path to the GenBank file.
    
    Returns:
    - organism_name (str): Name of the organism from the GenBank file.
    - section_name (str): Second-to-last term from the organism lineage.
    - gene_info (dict): Dictionary of gene names (in lowercase) associated with 'CDS' features and their sequences.
    """
    gene_info = {}
    organism_name = None
    section_name = None
    
    # Open and parse the GenBank file
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            organism_name = record.annotations.get('organism', 'Unknown').replace(" ", "_")
            # Extract the second-to-last term from the organism lineage
            lineage = record.annotations.get('taxonomy', [])
            if len(lineage) >= 2:
                section_name = lineage[-2]
            for feature in record.features:
                # Check for 'CDS' features and extract gene names and sequences
                if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0].lower()
                    sequence = feature.extract(record.seq)
                    gene_info[gene_name] = sequence
    
    return organism_name, section_name, gene_info

def process_genbank_files_for_genes(file_paths):
    """
    Process multiple GenBank files to extract gene names and sequences associated with 'CDS' features.
    
    Args:
    - file_paths (list): List of paths to GenBank files.
    
    Returns:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - all_genes (set): Set of all unique gene names found across all files.
    """
    organism_section_gene_map = {}
    all_genes = set()
    
    for file_path in file_paths:
        organism, section, genes = extract_genes_with_cds(file_path)
        organism_section_gene_map[organism] = {'section': section, 'genes': genes}
        all_genes.update(genes.keys())
    
    return organism_section_gene_map, all_genes

def write_gene_sequences(organism_section_gene_map, output_dir, overwrite):
    """
    Write the gene sequences to individual FASTA files in the output directory.
    
    Args:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - output_dir (str): Directory to write the gene sequence FASTA files.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    if os.path.exists(output_dir) and overwrite:
        shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    gene_sequences = {}
    
    for organism, info in organism_section_gene_map.items():
        for gene, sequence in info['genes'].items():
            if gene not in gene_sequences:
                gene_sequences[gene] = []
            gene_sequences[gene].append(SeqRecord.SeqRecord(sequence, id=organism, description=""))
    
    for gene, sequences in gene_sequences.items():
        gene_file = os.path.join(output_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file) or overwrite:
            SeqIO.write(sequences, gene_file, "fasta")
            print(f"Gene sequences for {gene} written to {gene_file}")

def align_gene_sequences(gene_sequences_dir, output_dir, overwrite):
    """
    Align the gene sequences using MAFFT and save the aligned sequences in the output directory.
    
    Args:
    - gene_sequences_dir (str): Directory containing gene sequence FASTA files.
    - output_dir (str): Directory to write the aligned gene sequence FASTA files.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    if os.path.exists(output_dir) and overwrite:
        shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for file_name in os.listdir(gene_sequences_dir):
        if file_name.endswith('.fasta'):
            input_path = os.path.join(gene_sequences_dir, file_name)
            output_path = os.path.join(output_dir, file_name)
            if not os.path.exists(output_path) or overwrite:
                # Run MAFFT alignment
                subprocess.run(f'mafft --threads -1 --auto {input_path} > {output_path}', shell=True)
                print(f"Aligned sequences for {file_name} written to {output_path}")

def run_raxml_ng_for_genes(aligned_sequences_dir, output_dir, overwrite):
    """
    Run RAxML-NG to generate phylogenetic trees with the GTR+GAMMA model and 100 bootstraps.
    
    Args:
    - aligned_sequences_dir (str): Directory containing aligned gene sequence FASTA files.
    - output_dir (str): Directory to write the RAxML-NG output files.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    if os.path.exists(output_dir) and overwrite:
        shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for file_name in os.listdir(aligned_sequences_dir):
        if file_name.endswith('.fasta'):
            input_path = os.path.join(aligned_sequences_dir, file_name)
            output_prefix = os.path.join(output_dir, file_name.replace('.fasta', ''))
            if not os.path.exists(f'{output_prefix}.raxml.bestTreeCollapsed') or overwrite:
                # Run RAxML-NG
                subprocess.run(f'raxml-ng --all --msa {input_path} --model GTR+G --bs-trees 100 --prefix {output_prefix}', shell=True)
                print(f"RAxML-NG analysis for {file_name} completed")

def perform_astral_analysis(raxml_trees_dir, output_dir, overwrite):
    """
    Run ASTRAL to generate a species tree from RAxML-NG best tree files with 100 multispecies bootstraps.
    
    Args:
    - raxml_trees_dir (str): Directory containing RAxML-NG best tree files.
    - output_dir (str): Directory to write the ASTRAL output files.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    
    # Ensure the output directory exists or remove and recreate it if overwrite is True
    if os.path.exists(output_dir) and overwrite:
        shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Combine RAxML best trees into a single file if not existing or overwrite flag is set
    best_tree_files = [os.path.join(raxml_trees_dir, f) for f in os.listdir(raxml_trees_dir) if f.endswith('.raxml.bestTreeCollapsed')]
    combined_tree_file = os.path.join(output_dir, 'in.tree')
    
    if not os.path.exists(combined_tree_file) or overwrite:
        with open(combined_tree_file, 'w') as outfile:
            for tree_file in best_tree_files:
                with open(tree_file, 'r') as infile:
                    outfile.write(infile.read())

    # Combine RAxML bootstrap files into a single list if not existing or overwrite flag is set
    bootstrap_files = [os.path.join(raxml_trees_dir, f) for f in os.listdir(raxml_trees_dir) if f.endswith('.raxml.bootstraps')]
    bootstrap_file_list = os.path.join(output_dir, 'bs-files')
    
    if not os.path.exists(bootstrap_file_list) or overwrite:
        with open(bootstrap_file_list, 'w') as outfile:
            for bs_file in bootstrap_files:
                outfile.write(bs_file + '\n')

    # Run ASTRAL with all support metrics (-t 2)
    astral_output_prefix = os.path.join(output_dir, 'astral_species_tree')
    if not os.path.exists(f'{astral_output_prefix}.tre') or overwrite:
        subprocess.run(f'astral -i {combined_tree_file} -b {bootstrap_file_list} -o {astral_output_prefix}.tre -t 2', shell=True)
        print(f"ASTRAL analysis completed. Output written to {astral_output_prefix}.tre")


def reorder_columns(df, section_order):
    """
    Reorder the columns of a DataFrame based on a user-specified list of section names.
    
    Args:
    - df (pd.DataFrame): DataFrame containing the gene summary.
    - section_order (list): List of section names in the desired order.
    
    Returns:
    - pd.DataFrame: DataFrame with reordered columns.
    """
    columns = df.columns.tolist()
    new_order = ['Gene']
    for section in section_order:
        new_order.extend([col for col in columns if section in col])
    # Add any remaining columns that were not in the section order
    new_order.extend([col for col in columns if col not in new_order])
    return df[new_order]

def write_summary_files(organism_section_gene_map, all_genes, output_base, feature_section_summary, section_order, overwrite):
    """
    Write the gene summary and section-wise summary to CSV and XLSX files.
    
    Args:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - all_genes (set): Set of all unique gene names found across all files.
    - output_base (str): Base path for the output files (without extension).
    - feature_section_summary (bool): Whether to generate the section-wise summary.
    - section_order (list): List of section names in the desired order.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    data = {'Gene': sorted(all_genes)}
    for organism in organism_section_gene_map:
        data[organism] = [gene if gene in organism_section_gene_map[organism]['genes'] else '' for gene in sorted(all_genes)]
    
    df = pd.DataFrame(data)
    
    # Reorder columns based on section order
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
    
    # Generate section-wise summary if requested
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
        
        # Reorder columns based on section order
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

def main():
    """
    Main function to parse command-line arguments and generate the gene summary files.
    """
    parser = argparse.ArgumentParser(description='Process GenBank files and extract gene names and sequences.')
    parser.add_argument('--input', nargs='+', required=True, help='Path to the GenBank files')
    parser.add_argument('--feature_section_summary', action='store_true', help='Generate section-wise feature summary')
    parser.add_argument('--gene_sequences', action='store_true', help='Generate gene sequence FASTA files')
    parser.add_argument('--align_sequences', action='store_true', help='Align gene sequences using MAFFT')
    parser.add_argument('--run_raxml', action='store_true', help='Run RAxML-NG for gene-tree calculations')
    parser.add_argument('--run_astral', action='store_true', help='Run Astral for coalescent-based phylogenetic tree inference')
    parser.add_argument('--output_dir', required=True, help='Base path for the output files (without extension)')
    parser.add_argument('--group_order', nargs='+', help='Order of sections for columns in the output')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files and directories')

    args = parser.parse_args()
    
    file_paths = args.input
    output_base = args.output_dir
    feature_section_summary = args.feature_section_summary
    gene_sequences = args.gene_sequences
    align_sequences = args.align_sequences
    run_raxml = args.run_raxml
    run_astral = args.run_astral
    section_order = args.group_order
    overwrite = args.overwrite
    
    # Process the GenBank files to extract gene information
    organism_section_gene_map, all_genes = process_genbank_files_for_genes(file_paths)
    
    # Write the summary to CSV and XLSX files
    write_summary_files(organism_section_gene_map, all_genes, output_base, feature_section_summary, section_order, overwrite)
    
    # Write gene sequences to FASTA files if requested
    if gene_sequences:
        gene_sequences_dir = os.path.join(output_base + "_gene_sequences")
        if not os.path.exists(gene_sequences_dir) or overwrite:
            write_gene_sequences(organism_section_gene_map, gene_sequences_dir, overwrite)
    
    # Align gene sequences using MAFFT if requested
    if align_sequences:
        aligned_sequences_dir = os.path.join(output_base + "_aligned_sequences")
        if not os.path.exists(aligned_sequences_dir) or overwrite:
            align_gene_sequences(gene_sequences_dir, aligned_sequences_dir, overwrite)
    
    # Run RAxML-NG analysis if requested
    if run_raxml:
        raxml_output_dir = os.path.join(output_base + "_raxml_trees")
        if not os.path.exists(raxml_output_dir) or overwrite:
            run_raxml_ng_for_genes(aligned_sequences_dir, raxml_output_dir, overwrite)
    
    # Run ASTRAL analysis if requested
    if run_astral:
        astral_output_dir = os.path.join(output_base + "_astral_trees")
        if not os.path.exists(astral_output_dir) or overwrite:
            perform_astral_analysis(raxml_output_dir, astral_output_dir, overwrite)

if __name__ == "__main__":
    main()
