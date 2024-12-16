# Tools-For-Finding-Comparing-Assembling-NCBI-Genomes
Automatically download, compare, and assemble genomes from NCBI


### - findGenome.py (version 5)

**Dependencies**
This script requires Python 3.6+ and the following libraries:
pip install biopython

Install the NCBI Datasets CLI for handling nuclear genomes:
conda install -c bioconda ncbi-datasets-cli

**This script will automatically download and filter organellar and nuclear genomes from the NCBI database**

```
# basic code:
python findGenome.py "group" "outfolder" --genome_type "chloroplast" --length_threshold INT --batch_size 50 --duplicate_removal --max_individuals_per_species INT --overwrite

# examples:
python findGenome5.py -g "ranunculaceae" -o ./genome_ranunculaceae --genome_type "nuclear_genome" --overwrite --email XXX@XXX

python findGenome5.py -g "ranunculaceae" -o ./chloroplast_ranunculaceae --genome_type "chloroplast" --duplicate_removal --max_individuals 2 --overwrite --email XXX@XXX

python findGenome5.py -g "ranunculaceae" -o ./mitogenome_ranunculaceae --genome_type "mitochondrial" --duplicate_removal --max_individuals 2 --overwrite --email XXX@XXX

# usage:
findGenome.py [-h] [--outfolfder OUTFOLDER] [--group GROUP]
                   [--genome_type {chloroplast,mitochondrial,nuclear_genome}] [--batch_size BATCH_SIZE]
                   [--duplicate_removal] [--max_individuals MAX_INDIVIDUALS_PER_SPECIES]
                   [--overwrite] [--email EMAIL]

Download plastid, mitochondrial, or nuclear genomes from NCBI.

options:
  -h, --help            Show this help message and exit
  -o, --outfolder       Output folder for downloaded files.
  -g, --group           Taxonomic group or organism name (e.g., genus, family, order).
  -t, --genome_type {chloroplast,mitochondrial,nuclear_genome}
                        The type of genome to download.
  --batch_size BATCH_SIZE
                        Batch size for downloading genomes (only organellar genomes).
  --duplicate_removal   Remove duplicate sequence files (only for organellar genomes).
  --max_individuals MAX_INDIVIDUALS_PER_SPECIES
                        Maximum number of individuals per species to retain (only organellar genomes).
  --overwrite           Overwrite existing output folder.
  --email               Your email for NCBI Entrez queries.
```

### - assembleGenome.py

**Dependencies**
This script requires Python 3.6+ and the following libraries:
pip install biopython

**This script automatically extracts and compares features, and aligns and assembles annotated CDS regions from sequence files downloaded from NCBI**

```
# basic code:
python assembleGenome.py -i INPUT -f FEATURE_SUMMARY -g --group_order GROUP_ORDER -s -a -r -t --output_dir OUTPUT_DIR --select_group SELECT_GROUP --overwrite

# example:
python assembleGenome.py -i ./plastomes_ranunculus/*.gb -f mitogenome_features -g -o Fumarioideae Thalictroideae Delphinieae Ranunculeae Anemoneae -s -a -r -t --output_dir ./phylogenies_ranunclaceae --overwrite

# usage:
assembleGenome.py [-h] -i INPUT [FILE1.gb FILE2.gb ...] [-f FEATURE_SUMMARY] [-g] [-o GROUP_ORDER [GROUP1 GROUP2 ...]]
                         [-s] [-a] [-r] [-t]
                         [--output_dir OUTPUT_DIR]
                         [--select_group SELECT_GROUP]
                         [--overwrite]

Process GenBank files and extract gene names and sequences.

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Paths to the GenBank files
  -f FEATURE_SUMMARY, --feature_summary FEATURE_SUMMARY
                        Base path for feature summary files (CSV and XLSX)
  -g, --group_feature_summary
                        Whether to generate a section-wise feature summary
  -o GROUP_ORDER [GROUP_ORDER ...], --group_order GROUP_ORDER [GROUP_ORDER ...]
                        Optional: Order of sections in the summary files
  -s, --generate_gene_sequences
                        Generate gene sequences in FASTA format
  -a, --align_sequences
                        Align gene sequences using MAFFT
  -r, --run_raxml       Run RAxML-NG for phylogenetic analysis
  -t, --run_astral      Perform ASTRAL analysis on RAxML trees
  --select_group SELECT_GROUP
                        Limit gene extraction to a specific group
  --output_dir OUTPUT_DIR
                        Output directory for gene sequences and alignment results
  --overwrite           Overwrite existing files if they exist
```

### If you use any of the scripts, please cite the following reference until the journal article is published: 
Karbstein et al. (2023), BioRxiv (https://doi.org/10.1101/2023.08.08.552429)
