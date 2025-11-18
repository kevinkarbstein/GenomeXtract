# Find, Compare, and Assembly NCBI-Genomes
Automatically download, compare, and assemble genomes from NCBI

**Dependencies**
This script requires Python 3.6+ and the following libraries:

conda install python=3.10 

pip install biopython

pip install pandas

Install the NCBI Datasets CLI for handling nuclear genomes:

conda install -c bioconda ncbi-datasets-cli

conda install bioconda::raxml-ng # v1.2.2

conda install bioconda::iqtree # v3.0.1

conda install bioconda::mafft # v7.525


### - findGenome.py

**This script will automatically download and filter organellar and nuclear genomes from the NCBI database**

```
# basic code:
python findGenome.py --group "group" --outfolder "outfolder" --genome_type "chloroplast" --batch_size 50 --duplicate_removal --max_individuals_per_species INT --overwrite

# examples:
python findGenome.py --group "ranunculaceae" --outfolder ./genome_ranunculaceae --genome_type "nuclear_genome" --overwrite --email XXX@XXX

python findGenome.py --group "ranunculus" --outfolder ./chloroplast_ranunculus --genome_type "chloroplast" --duplicate_removal --max_individuals 2 --batch_size 50 --overwrite --email XXX@XXX

python findGenome.py --group "ranunculaceae" --outfolder ./mitogenome_ranunculaceae --genome_type "mitochondrial" --duplicate_removal --max_individuals 2 --batch_size 50 --overwrite --email XXX@XXX

# usage:
findGenome.py [-h] [--outfolfder OUTFOLDER]
                   [--group GROUP]
                   [--genome_type {chloroplast,mitochondrial,nuclear_genome}]
                   [--annotated]
                   [--assembly_level {scaffold, chromosome}]
                   [--batch_size BATCH_SIZE]
                   [--duplicate_removal]
                   [--max_individuals MAX_INDIVIDUALS]
                   [--overwrite]
                   [--email EMAIL]

Download plastid, mitochondrial, or nuclear genomes from NCBI.

options:
  -h, --help            Show this help message and exit
  --outfolder           Output folder for downloaded files (STRING).
  --group               Taxonomic group or organism name (STRING; e.g., the name of the genus, family, order).
  --genome_type         The type of genome to download (STRING; chloroplast,mitochondrial,nuclear_genome).
  --annotated           Select only gene-annotated nuclear genomes.
  --assembly_level      Choose the assmbly level of the nuclear genome (STRING; scaffold, chromosome).
  --batch_size          Batch size for downloading genomes (INT; only organellar genomes).
  --duplicate_removal   Remove duplicate sequence files (only for organellar genomes). Prioritize NC_* or the latest release(s).
  --max_individuals     Maximum number of individuals per species to retain (INT; only organellar genomes). Prioritize NC_* the latest release(s).
  --overwrite           Overwrite existing output folder.
  --email               Your email for NCBI Entrez queries (STRING).
```

### - findClosestGenome.py

**This script will automatically download and filter organellar and nuclear genomes from the NCBI database**

```
# basic code:
python findClosestGenome.py --taxon "species/taxon" "outfolder" --genome_type "chloroplast" --overwrite

# examples:
python findClosestGenome.py --taxon "ranunculus cassubicifolius" --outfolder ./closest_genome_ranunculaceae --genome_type "nuclear_genome" --overwrite --email XXX@XXX

python findClosestGenome.py --taxon "ranunculus" --outfolder ./closest_chloroplast_ranunculus --genome_type "chloroplast" --overwrite --email XXX@XXX

python findClosestGenome.py --taxon "ranunculus" --outfolder ./closest_mitogenome_ranunculaceae --genome_type "mitochondrial" --overwrite --email XXX@XXX

# usage:
findGenome.py [-h] [--outfolfder OUTFOLDER]
                   [--taxon SPECIES]
                   [--genome_type {chloroplast,mitochondrial,nuclear_genome}]
                   [--annotated]
                   [--assembly_level {scaffold, chromosome}]
                   [--overwrite]
                   [--email EMAIL]

Find the closest available reference genomes(s) of a given taxon in NCBI.

options:
  -h, --help            Show this help message and exit
  --outfolder           Output folder for the result file (STRING).
  --taxon               Species or higher-level taxon name (STRING; e.g., Genus or Family).
  --genome_type         The type of genome to download (STRING; chloroplast,mitochondrial,nuclear_genome).
  --annotated           Select only gene-annotated nuclear genomes.
  --assembly_level      Choose the assmbly level of the nuclear genome (STRING; scaffold, chromosome).
  --overwrite           Overwrite existing output folder.
  --email               Your email for NCBI Entrez queries (STRING).
```

### - assembleGenome.py / assembleGenes.py

**These scripts automatically extract and compare features, and align and assemble sequences/annotated CDS regions from files downloaded from NCBI**

```
# basic code:
python assembleGenome.py -i INPUT --feature_summary FEATURE_SUMMARY --feature_summary --group_order GROUP_ORDER --generate_gene_sequences --align_sequences --run_raxml --output_dir OUTPUT_DIR --select_group SELECT_GROUP --overwrite

# example:
For a small test datasets, run: python findGenome.py --group "ranunculus" --outfolder ./chloroplast_ranunculus --genome_type "chloroplast" --duplicate_removal --max_individuals 2 --overwrite --email XXX@XXX
python assembleOrgGenome.py -i chloroplast_ranunculus/*.gb --feature_summary ranunculus_features --generate_gene_sequences --align_sequences --run_raxml --output_dir chloroplast_ranunculus_raxml_ng/
python assembleOrgGenome.py -i chloroplast_ranunculus/*.gb --feature_summary ranunculus_features --generate_gene_sequences --align_sequences --run_iqtree --output_dir chloroplast_ranunculus_iqtree/

# usage:
assembleOrgGenome.py [-h] [--input INPUT [FILE1.gb FILE2.gb ...]
                       [--feature_summary FEATURE_SUMMARY]
                       [--group_feature_summary]
                       [--group_order GROUP_ORDER [GROUP1 GROUP2 ...]]
                       [--generate_gene_sequences]
                       [--align_sequences]
                       [--run_raxml]
                       [--raxml_model RAXML_MODEL]
                       [--bootstrap_replicates BOOTSTRAP_REPLICATES]
                       [--run_iqtree]
                       [--output_dir OUTPUT_DIR]
                       [--select_group SELECT_GROUP]
                       [--overwrite]

Process GenBank files and extract gene names and sequences.

options:
  -h, --help                show this help message and exit
  --input INPUT             Path to the GenBank files (STRING)
  --feature_summary         Folder name for feature summary files (STRING)                  
  --group_feature_summary   Whether to generate a section-wise feature summary               
  --group_order             Order of sections in the summary files (STRING)          
  --generate_gene_sequences Generate gene sequences in FASTA format       
  --align_sequences         Align gene sequences using MAFFT
  --run_raxml               Run RAxML-NG for phylogenetic analysis
  --raxml_model             RAxML-NG model to use (STRING; default: GTR+G)
  --bootstrap_replicates    Number of bootstrap replicates for RAxML-NG (INT; default: 100)
  --run_iqtree              Run IQ-TREE for phylogenetic analysis
  --select_group            Limit gene extraction to a specific group (STRING)  
  --output_dir              Output directory for gene sequences and alignment results (STRING)
  --overwrite               Overwrite existing files if they exist
```



```
# basic code:
python assembleGenes.py -i INPUT -s FEATURE_SUMMARY -g --group_order GROUP_ORDER -s -a -r --output_dir OUTPUT_DIR --select_group SELECT_GROUP --overwrite

# example:
For a small test datasets, run: python findGenome5.py -g "ranunculus" -o ./chloroplast_ranunculus --genome_type "chloroplast" --duplicate_removal --max_individuals 2 --overwrite --email XXX@XXX
python assembleGenes.py -i chloroplast_ranunculus/*.gb -o ranunculus_gene_features -s -g -a -r -x 

For a small test datasets, run: python findGenome5.py -g "ranunculales" -o ./mitogenome_ranunculales --genome_type "mitochondrial" --duplicate_removal --max_individuals 2 --overwrite --email XXX@XXX
python assembleGenes.py -i mitogenome_ranunculales/*.gb -o ranunculales_gene_features —feeature_section_summary -o Papaveroideae Fumarioideae Thalictroideae Delphinieae Ranunculeae Anemoneae -s -g -a -r -x

# usage:
assembleGenes.py [-h] --input INPUT [FILE1.gb FILE2.gb ...]
                 [-o GROUP_ORDER [GROUP1 GROUP2 ...]]
                 [--feature_section_summary]
                 [--generate_gene_sequences]
                 [--align_sequences]
                 [--run_raxml]
                 [--run_astral]
                 [--output_dir OUTPUT_DIR]
                 [--select_group SELECT_GROUP]
                 [--overwrite]

Process GenBank files and extract gene names and sequences.

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Path to the GenBank files
  -f, --feature_section_summary
                        Generate section-wise feature summary
  --group_order GROUP_ORDER [GROUP_ORDER ...]
                        Limit gene extraction to a specific group
  -g, --generate_gene_sequences
                        Generate gene sequences in FASTA format
  -a, --align_sequences
                        Align gene sequences using MAFFT
  -r, --run_raxml       Run RAxML-NG for phylogenetic analysis
  -o, --output_dir OUTPUT_DIR
                        Output directory for gene sequences and alignment results
  --overwrite           Overwrite existing files if they exist
```


### If you use any of the scripts, please cite the following article: 
**Karbstein K, Choudhary N, Xie T, Tomasello S, Wagner ND, Barke BH, Paetzold C, Bradican JP, Preick M, Himmelbach A, Stein N, Papantonis A, Irisarri I, de Vries J, Pucker B, Hörandl E**. Assembling genomes of non-model plants: A case study with evolutionary insights from Ranunculus (Ranunculaceae). **The Plant Journal**, 2025; 00; e70390. https://doi.org/10.1111/tpj.70390
