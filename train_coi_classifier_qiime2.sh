#!/usr/local/bin
# The code snippet to train the bayesian classifier
# We will have the sequences, taxonomy and classifier in the same folder

# usage bash train_classifier_qiime2.sh <path_to_folder>

FOLDER=$1

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${FOLDER}"/midori_sequences.qza \
  --i-reference-taxonomy "${FOLDER}"/midori_taxonomy.qza \
  --o-classifier "${FOLDER}"/midori-coi-classifier.qza
  
  
# Inspect the classifier

qiime tools peek "${FOLDER}"/midori-coi-classifier.qza