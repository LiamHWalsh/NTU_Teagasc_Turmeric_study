#!/bin/sh
#SBATCH --job-name 01_fasnicar_Nottingham
#SBATCH --error 01_fasnicar_Nottingham.err
#SBATCH --output 01_fasnicar_Nottingham.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=liam.walsh@teagasc.ie 

# Define input, output, and index directory.
#in_dir="/data/Food/analysis/R1838_DOMINO/MASTER/00_download"
out_dir="/data/Food/analysis/R1909_CODIET/NottinghamJuly2025-460410994/01_preprocessing/01_fasnicar"
index_dir="/data/Food/analysis/R1838_DOMINO/MASTER/01_preprocessing_v2/01_bowtie_index"



folder="/data/Food/primary/R1909_CODIET/NottinghamJuly2025-460410994/"

rm -r $out_dir

        # Print only the last field after the last delimiter "/" in the path and select only \
#the base name from the in_dir.
        folder_base=$(echo "$folder" | awk -F / '{print $NF}')

 folder_base=$(echo "$folder" | awk -F / '{print $NF}')

# Find all the files *.fastq.gz from the in_dir; Print only the last field after the \
#last delimiter "/"
        # Remove anything that contains _..fastq.gz and other formats with nothing that will \
#leave only the base name; Avoid duplication.
  
for filename in $(find "$folder" -type f -name "*.fastq.gz" | awk -F / '{print $NF}' | sed 's/_..fastq.gz//g' | sed 's/_L00.*//g' | sort -u)
do
        # Modify the filename by removing the SRR or ERR in the base name.
        filename_modify=$(echo $filename | sed 's/SRR//;s/ERR//')

        # Make an output directory that contains the modified file name.
        mkdir -p $out_dir/$folder_base/$filename_modify

        # Copy the modified directory to the output directory.
         # error is here 
          find "$folder" -type f \( -name "*$filename*" -and -name "*.fastq.gz" \)  -exec cp {} $out_dir/$folder_base/$filename_modify \; 

# Load the preprocessing tool.
module load fasnicar/0.2.4

        # Run preprocessing script.
        preprocess.new.py -e .fastq.gz -f _R1 -r _R2 -x "$index_dir" --verbose -n "$SLURM_CPUS_\
PER_TASK" -i "$out_dir"/"$folder_base"/"$filename_modify"

# Unload the preprocessing tool.
module unload fasnicar/0.2.4

done

done

