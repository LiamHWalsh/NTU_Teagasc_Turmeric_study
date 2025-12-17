#!/bin/sh
#SBATCH --error err_metagenomics_pipeline
#SBATCH --output out_metagenomics_pipeline
#SBATCH --job-name metagenomics_pipeline
#SBATCH -p Priority,Background,GPU
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=liam.walsh@teagasc.ie

# Load required modules
module load nextflow/24.04.2_5914
module load apptainer/1.1.9

# Run the metagenomics pipeline with HUMAnN3
nextflow run main_metagenomics_v3.nf \
  --input_dir /data/Food/primary/R1909_CODIET/NottinghamJuly2025-460410994/ \
  --output_dir /data/Food/analysis/R1909_CODIET/NottinghamJuly2025-460410994/results_v2/ \
  --host_index_dir /data/Food/analysis/R1838_DOMINO/MASTER/01_preprocessing_v2/01_bowtie_index \
  --metaphlan_db_dir /data/Food/analysis/R1838_DOMINO/MASTER/04_short_read_taxonomic_profiling/04metaphlan/04metaphlan_db/ \
  --metaphlan_index mpa_vOct22_CHOCOPhlAnSGB_202403 \
   --humann_nucleotide_db /data/Food/analysis/R1909_CODIET/NottinghamJuly2025-460410994/databases/humann_db/chocophlan \
  --humann_protein_db /data/Food/analysis/R1909_CODIET/NottinghamJuly2025-460410994/databases/humann_db/uniref \
  --cpus 8 \
  --email liam.walsh@teagasc.ie \
  -resume \
  -process.executor slurm


# Unload modules
module unload nextflow/24.04.2_5914
module unload apptainer/1.1.9

echo "Pipeline submission completed"
