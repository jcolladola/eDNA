#!/bin/bash
#SBATCH --job-name=taxonomy_array
#SBATCH --output=/home/meg/rgallego/reports/taxonomy_%A_%a.out
#SBATCH --error=/home/meg/rgallego/reports/taxonomy_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00


## BEFORE RUNNING THIS JOB:
## You need to use tazonkit to retrive the lineages of all taxids
##  awk '{print $13}' blast_output_BOLD/split_part_00*  | sort | uniq | taxonkit  reformat -I 1 -f "{k};{p};{c};{o};{f};{g};{s}" > lineage.txt
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
Rscript "${SCRIPT_DIR}"/parse_blast_taxonkit.r "$BLAST_FILE"

