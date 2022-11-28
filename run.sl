#!/bin/bash -e
#SBATCH --job-name=crossNational  # job name (shows up in queue)
#SBATCH --time=01-00:00:00        # Walltime (DD-HH:MM:SS)
#SBATCH --mem=50000               # total memory in MB
#SBATCH --cpus-per-task=4         # 4 CPUs
#SBATCH --account=uoa03415        # Project code

# load R
module load R/4.2.1-gimkl-2022a

# load necessary packages
module load MPFR/4.1.0-GCC-11.3.0
module load UDUNITS/2.2.26-GCCcore-11.3.0
module load GCCcore/5.4.0
module load GDAL/3.5.1-gimpi-2022a

# run script
Rscript run.R
