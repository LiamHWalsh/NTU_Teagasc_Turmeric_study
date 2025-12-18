#!/bin/sh
#SBATCH --error err_test
#SBATCH --output out_test
#SBATCH --job-name test
#SBATCH -p Priority,Background,GPU
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=liam.walsh@teagasc.ie

# ==== Load Apptainer ====
module load apptainer/1.1.9

# ==== Define directories ====
R_dir=/data/Food/analysis/R1838_DOMINO/MASTER/containers
WORKDIR=/data/Food/analysis/R1909_CODIET/liam_walsh/r

cd $WORKDIR || exit 1

# ==== Run R inside container ====
apptainer exec \
  --bind $R_dir/Rlibs:/usr/local/lib/R/site-library \
  --bind /data:/data \
  --bind /usr/bin/git:/usr/bin/git \
  $R_dir/bioconductor_docker_RELEASE_3_19.sif \
  Rscript differnetial_abundance_notthingham.R

