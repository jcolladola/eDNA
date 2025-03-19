#!/usr/local/bin/bash

# USAGE   bash vsearch.sh <folder>
SCRIPT_folder=$(dirname $0)
INPUTfolder=$1





echo "INPUT_folder; ${INPUTfolder}"


vsearch  --uchime_denovo "${INPUTfolder}"/vsearch_input.fasta --nonchimeras "${INPUTfolder}"/non_chim.fasta --chimeras "${INPUTfolder}"/chim.fasta
