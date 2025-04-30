#!/bin/bash
#SBATCH --job-name=taxonomy_array
#SBATCH --output=/home/meg/rgallego/reports/taxonomy_%A_%a.out
#SBATCH --error=/home/meg/rgallego/reports/taxonomy_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

## TO LAUNCH THIS JOB:
## out_folder=/home/meg/rgallego/output
## lines=$(ls "$out_folder"/split_part* | wc -l)

## sbatch --array=1-${lines} --export=OUT_FOLDER="$out_folder" taxonomy_array.sh


# Load modules
# module load R # so far I do not use modules 
SCRIPT_DIR="/home/meg/rgallego/ccc_cluster_scripts/scripts/taxonomy"
# Get the folder passed via sbatch
IFS=$'\n' BLAST_FILES=($(ls -1v "${OUT_FOLDER}"/split_part*))
BLAST_FILE=${BLAST_FILES[$SLURM_ARRAY_TASK_ID - 1]}

echo "Processing $BLAST_FILE"
Rscript "${SCRIPT_DIR}"/parse_blast_rCRUX.r "$BLAST_FILE"

## After all is done, you can combine all outputs in R with
# list.files("your_folder", pattern = "max_pidents_.*\\.tsv", full.names = TRUE) |>
#  map_dfr(read_tsv) -> all_max_pidents
