#!/usr/local/bin/bash

# USAGE   bash vsearch.sh <folder with centroids output> <path/to/vsearchoutput>
SCRIPT_folder=$(dirname $0)
INPUTfolder=$1
OUTPUTfolder="${2}"

mkdir "${OUTPUTfolder}"
mkdir "${OUTPUTfolder}"/global
mkdir "${OUTPUTfolder}"/bysample


echo "INPUT_folder; ${INPUTfolder}"

shopt -s nullglob
fasta_files=("${INPUT_folder}"/*.fasta)
shopt -u nullglob


echo "Found files: ${fasta_files[@]}"


# for fasta_file in "${INPUTfolder}"/*.fasta; do

#     sample=$(echo "${fasta_file%.*}")
#     output=$(basename "${sample}")

#     vsearch  --uchime_denovo "${fasta_file}" --nonchimeras "${OUTPUTfolder}"/bysample/"${output}"_non_chim.fasta --chimeras "${OUTPUTfolder}"/bysample/"${output}"_chim.fasta

# done

# Global search

# Rscript "${SCRIPT_folder}"/merge.fastas.R "${INPUTfolder}" "${OUTPUTfolder}"

vsearch  --uchime_denovo "${OUTPUTfolder}"/all_together.fasta --nonchimeras "${OUTPUTfolder}"/global/non_chim.fasta --chimeras "${OUTPUTfolder}"/global/chim.fasta
