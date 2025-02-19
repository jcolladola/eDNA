#!/usr/local/bin

# USAGE   bash swarm.sh <folder with all demultiplexed fastas>

folder=$1
CENTROIDS="${folder}"/centroids
OUTPUTS="${folder}"/outputs
mkdir "${CENTROIDS}"
mkdir "${OUTPUTS}"

for fasta_file in "${folder}"/*fasta; do

    sample="${fasta_file%.*}"
    output=$(basename "${sample}")

    swarm -d 2 -t 4 -z -w "${CENTROIDS}"/"${output}".centroids.fasta "${fasta_file}"   -o "${OUTPUTS}"/"${output}".output.txt

done