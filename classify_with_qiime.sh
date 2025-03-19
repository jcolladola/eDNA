#!/usr/local/bin
# The code snippet to classify a fasta file with the QIIME bayesian classifier
# We will have the sequences and classifier as arguments

# usage bash train_classifier_qiime2.sh <path_to_fasta> <path_to_classifier> <OUTPUT_FOLDER>


INPUT_FASTA=$1

CLASSIFIER=$2

OUTPUT_FOLDER=$3

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "${INPUT_FASTA}" \
  --output-path "${OUTPUT_FOLDER}"/sequences.qza


# Asigna taxonomía a las secuencias representativas utilizando el clasificador entrenado
qiime feature-classifier classify-sklearn \
--i-classifier "${CLASSIFIER}" \
--i-reads "${OUTPUT_FOLDER}"/sequences.qza \
--o-classification "${OUTPUT_FOLDER}"/classified_seqs.qza \
--p-n-jobs 1 \
--p-reads-per-batch 10000 \
--p-pre-dispatch 1 \
--verbose 