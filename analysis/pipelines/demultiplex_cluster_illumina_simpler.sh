#!/usr/bin/env bash

# Usage bash demultiplex_both_fastqs.sh banzai_params.sh
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template
 module load cutadapt/4.1
 module load miniconda/3.11
# Weird
#We need to gather: Location of functions  and fastqs:
MAIN_DIR="$(dirname "$0")"
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/* ; do
	source "${file}"
done

param_file=${1}

echo "Reading analysis parameters from:"
echo "${param_file}"
source "${param_file}"

# Check if the metadata file exists
if [[ -s "${SEQUENCING_METADATA}" ]]; then
	echo "Reading metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo 'ERROR! Could not find metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# Now fix line ends if needed

if [[ $( file "${SEQUENCING_METADATA}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${SEQUENCING_METADATA%.*}"

  EXT="${SEQUENCING_METADATA##*.}"

  NEWLINES_FIXED="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${SEQUENCING_METADATA}" > "${NEWLINES_FIXED}"

  echo "the old file was: ${SEQUENCING_METADATA}"

  echo "The new file is here:"

  echo "${NEWLINES_FIXED}"

else

  echo "The file passes test for CRLF. Everybody dance!"
  echo

fi

if [[ -s "${NEWLINES_FIXED}" ]]; then
	SEQUENCING_METADATA="${NEWLINES_FIXED}"
fi

OUTPUT_DIRECTORY_TEMP=/temporal/$SLURM_JOB_USER/$SLURM_JOB_ID

#Create output directory
START_TIME=$(date +%Y%m%d_%H%M)
OUTPUT_DIR="${OUTPUT_DIRECTORY_TEMP}"/demultiplexed_"${START_TIME}"
OUTPUT_DIR="${OUTPUT_DIRECTORY}"/demultiplexed_"${START_TIME}"
OUTPUT_DIRECTORY_TEMP="${OUTPUT_DIRECTORY}"
mkdir "${OUTPUT_DIR}"
echo "Output directory is ${OUTPUT_DIR}"
mkdir "${OUTPUT_DIR}"/input

# copy metadata and parameters file to output directory
cp "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/metadata.csv
cp "${param_file}" "${OUTPUT_DIR}"/banzai_params.sh

# Write a log file
LOGFILE="${OUTPUT_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1


DEMULT_DIR="${OUTPUT_DIR}"/demultiplexed
mkdir "${DEMULT_DIR}"
################################################################################
# READ METADATA
################################################################################
# report metadata dimensions
METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))
echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
## NOW WE HAVE LOADED THE SEQUENCING_METADATA - WE NEED to find the columns specified
## in the params file. We should set up an alert & quit if a critical column is not found

# Filnames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")
# Pass check
# Library names
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")

# COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

# Secondary indices
COLNUM_ID2=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")

# Secondary index sequence positions
# COLNUM_ID2_START=$( get_colnum "${COLNAME_ID2_START}" "${SEQUENCING_METADATA}")

# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")

# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_FILE2 COLNUM_ID1 COLNUM_ID2 \
COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2)

echo "Checking that all columns in metadata are there"

for column in "${all_columns[@]}" ; do

 if [ "${!column}" -gt 0 ]; then
	 echo "looking good, ${column}"
 else
  echo "Something went wrong with column name ${column}"
	echo "exiting script"
	exit
fi
done
echo "All columns passed test"



################################################################################
# CHECK FILES
################################################################################

if [[ "${ALREADY_DEMULTIPLEXED}" != "YES" ]]; then


	FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))

	FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
	  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
	  sort | uniq ))

	NFILE1="${#FILE1[@]}"
	NFILE2="${#FILE2[@]}"
	if [ "${NFILE1}" != "${NFILE2}" ]; then
		echo "ERROR: Whoa! different number of forward and reverse files"
	fi

	if [[ -n "${FILE1}" && -n "${FILE2}" ]]; then
	  echo 'Files read from metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
	  echo 'File names:'
		for (( i=0; i < "${NFILE1}"; ++i)); do
			printf '%s\t%s\n' "${FILE1[i]}" "${FILE2[i]}"
		done
		echo
	else
	  echo 'ERROR:' 'At least one file is not valid'
	  echo 'Looked in metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
	  echo 'Aborting script'
	  exit
	fi
	#here we play again
	if [[ "${SECONDARY_INDEX}" == "YES" ]]; then

		ID2S=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
		  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
		  sort | uniq))
		N_index_sequences="${#ID2S[@]}"
		ID2_LENGTH=${#ID2S[0]}
		# ID2_START=($(awk -F',' -v COLNUM=$COLNUM_ID2_START \
		#   'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
		#   sort | uniq))

		# check if number of indexes is greater than one:
		if [[ "${N_index_sequences}" -gt 1 ]]; then
			echo "Secondary indexes read from sequencing metadata (""${N_index_sequences}"" total)"
			echo
		else
		  echo
		  echo 'ERROR:' "${N_index_sequences}" 'index sequences found. There should probably be more than 1.'
		  echo
		  echo 'Aborting script.'
			exit
		fi

	fi
	echo "These are the secondary barcodes"
	echo "${ID2S[@]}"
	echo "that is, ${#ID2S[@]} unique barcodes"
	echo "and they seem to be sorted alphabetically?"
	echo "and they are this long "
	echo "ID2_LENGTH  es ${ID2_LENGTH}"


	################################################################################
	# Read in primers
	################################################################################
	PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq ))

	PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq ))

	if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
	  echo 'Primers read from metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
	  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
		echo
	else
	  echo 'ERROR:' 'At least one primer is not valid'
	  echo 'Looked in metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
	  echo 'Aborting script'
	  exit
	fi




	#######
	#Unique samples are given by combining the primary and secondary indexes
	######
	ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2 \
	'NR>1 {
	  print ";ID1=" $COLNUM1 ";ID2=" $COLNUM2
	}' "${SEQUENCING_METADATA}" )

	SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))

	#####
	# Check that sample names are not repeated
	#####
	NSAMPLES="${#SAMPLE_NAMES[@]}"

	# Now calculate the number of unique sample names 
	UNIQ_SAMPLES=( $(echo "${SAMPLE_NAMES[@]}" | tr ' ' '\n' | sort -u))
	N_UNIQ_SAMPLES="${#UNIQ_SAMPLES[@]}"


	if [[ "${NSAMPLES}" != "${N_UNIQ_SAMPLES}" ]]; then
		echo " At least one sample name is repeated "
		echo " I am not angry, just dissapointed. Exiting script"
		exit
	fi


	ID1_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
	ID1S=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}"  |\
			sort | uniq))
	ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
	ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))

	# write file for translating demultiplexed output to samples
	SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
	for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do
	  printf "ID1=%s;ID2A=%s;ID2B=%s\t%s_%s\t%s\n" \
		"${ID1_ALL[i]}" "${ID2_ALL[i]}" "${ID2_ALL_RC[i]}" \
		"${ID1_ALL[i]}" "${ID2_ALL[i]}" \
		"${SAMPLE_NAMES[i]}" >> "${SAMPLE_TRANS_FILE}"
	done
	for (( i=0; i < "${#ID1S[@]}"; i++ )); do
	  printf "File1:%s\tFile2:%s\tLib:%s\n" \
	  "${FILE1[i]}" "${FILE2[i]}" "${ID1S[i]}"


	done


	PRIMERS_FILE="${OUTPUT_DIR}"/pcr_primers.fasta
	PRIMERS_R2_FILE="${OUTPUT_DIR}"/pcr_primers_r2.fasta

	printf ">FWD\n${PRIMER1}\n>REV\n${PRIMER2}\n" > "${PRIMERS_FILE}"

	printf ">REV\n${PRIMER2}\n>FWD\n${PRIMER1}\n" > "${PRIMERS_R2_FILE}"

	source "${SCRIPT_DIR}"/functions/check_primers.sh "${PRIMERS_FILE}"

	source "${SCRIPT_DIR}"/functions/check_primers.sh "${PRIMERS_R2_FILE}"
	#Hooray it works
	#Create a dir for all the demultiplexed files

	# now we have to remove all hard-coded stuff and link it to
	#banzai_params
	# to get the .1 files trimmed and the .2 selected along
	# 
	# 	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	# 	printf "First_trim_R1,nReads_R1,First_trim_R2,nReads_R2,Second_trim_R1,nReads_StR1,Second_trim_R2,nReads_StR2,NEW_OUTPUT_Fwd_1,nreadsFwd,NEW_OUTPUT_Rev_1,nreadsRev\n" \
	# 	> "${OUTPUT_SUMMARY}"
	
	# long table version of the output summary
	
		OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
		printf "fastq_header,name,value,step,Sample\n" \
		> "${OUTPUT_SUMMARY}"
	

	################################################################################
	# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
	################################################################################

	for (( i=0; i < "${#FILE1[@]}"; i++ )); do
	  # Identify the forward and reverse fastq files.

	  # copy the two files to the temp file
		#   NEW_parent=/temporal/$SLURM_JOB_USER/$SLURM_JOB_ID/input
		#   mkdir "${NEW_parent}"
		
	  READ1="${PARENT_DIR}/${FILE1[i]}"
	  READ2="${PARENT_DIR}/${FILE2[i]}"

	  BASE1="${FILE1[i]%.*}"
	  BASE2="${FILE2[i]%.*}"

	  LIBNAME=$(awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v NAME=$COLNUM_ID1 \
	'{if ($COLNUM == VALUE) { print $NAME } }' $SEQUENCING_METADATA | sort | uniq)

	  mkdir "${OUTPUT_DIR}"/"${LIBNAME}"

	  INFO_FILE1="${OUTPUT_DIR}"/"${LIBNAME}"/info_"${LIBNAME}"_R1.txt
	  INFO_FILE2a="${OUTPUT_DIR}"/"${LIBNAME}"/info_"${LIBNAME}"_R2_a.txt
	  INFO_FILE2b="${OUTPUT_DIR}"/"${LIBNAME}"/info_"${LIBNAME}"_R2_b.txt
	  INFO_FILE2="${OUTPUT_DIR}"/"${LIBNAME}"/info_"${LIBNAME}"_R2.txt

	  DEMULT_FILE1="${OUTPUT_DIR}"/"${LIBNAME}"/info_demult_"${BASE1}".txt
	  DEMULT_FILE2="${OUTPUT_DIR}"/"${LIBNAME}"/info_demult_"${BASE2}".txt

	  echo "Info file"
	  echo "${INFO_FILE1}"

	## making two primers file for this library


	
	## making a barcode file for this library

	BARCODES_FILE="$OUTPUT_DIR"/barcodes_"${LIBNAME}".fasta
	awk -F',' -v COLNUM=$COLNUM_FILE1 -v VALUE=${FILE1[i]} -v ADAP=$COLNUM_ID2 \
	'{if ($COLNUM == VALUE) { printf ">%s\n%s\n", $ADAP, $ADAP } }' $SEQUENCING_METADATA > "${BARCODES_FILE}"

	BARCODES_FILE_R2="$OUTPUT_DIR"/barcodes_"${LIBNAME}"_R2.fasta

	cp "${BARCODES_FILE}" "${BARCODES_FILE_R2}" 

	  
	  mkdir "${OUTPUT_DIR}"/reports
	  mkdir "${OUTPUT_DIR}"/temp

		echo "Working on Library $[i+1] out of ${#FILE1[@]}"

	## First cutdapt: demultiplex both fastqs at once

	cutadapt -g "file:"${BARCODES_FILE}";min_overlap=8" -G "file:"${BARCODES_FILE_R2}";min_overlap=8"  --cores=0 --no-indels --pair-adapters \
	 "${READ1}" "${READ2}"\
	  -o "${OUTPUT_DIR}"/temp/"${LIBNAME}"_{name}.R1.fastq -p "${OUTPUT_DIR}"/temp/"${LIBNAME}"_{name}.R2.fastq \
	  -m 1 --json "${OUTPUT_DIR}"/reports/report_.1."${LIBNAME}".json --discard-untrimmed > "${OUTPUT_DIR}"/reports/"${LIBNAME}"_demult_logfile.txt


	inputReads=$(grep "Total read pairs processed" "${OUTPUT_DIR}"/reports/"${LIBNAME}"_demult_logfile.txt  | awk '{print $5}')
	echo "${BASE1},${LIBNAME},${inputReads},input,${LIBNAME}" >>  "${OUTPUT_SUMMARY}"
	

	## Loop through all demult_files, remove primers

	for demult_file in "${OUTPUT_DIR}"/temp/"${LIBNAME}"_*.R1.fastq; do

		demult_file2=$(echo ${demult_file} | sed "s/.R1.fastq/.R2.fastq/g" )

		base=$(basename "${demult_file}"  )
		base_no_ext="${base%.*.*}"
		barcode=$(echo "${base_no_ext}" | awk -F '_' '{print $NF}' )
		sample="${LIBNAME}"_"${barcode}"

		cutadapt -g "file:"${PRIMERS_FILE}";min_overlap=26" -G "file:"${PRIMERS_R2_FILE}";min_overlap=26" --cores=0 --no-indels --pair-adapters \
	 		"${demult_file}" "${demult_file2}" \
	  	-o "${DEMULT_DIR}"/"${sample}"_{name}.1.fastq -p "${DEMULT_DIR}"/"${sample}"_{name}.2.fastq \
	  	-m 1 --json "${OUTPUT_DIR}"/reports/report_.1."${LIBNAME}".json --discard-untrimmed > "${OUTPUT_DIR}"/reports/"${LIBNAME}"_primer_logfile.txt

		inputReads=$(grep "Total read pairs processed" "${OUTPUT_DIR}"/reports/"${LIBNAME}"_primer_logfile.txt  | awk '{print $5}')
		echo "${BASE1},${LIBNAME},${inputReads},demult,${sample}" >>  "${OUTPUT_SUMMARY}"

		outputReads=$(grep "Pairs written (passing filters):" "${OUTPUT_DIR}"/reports/"${LIBNAME}"_primer_logfile.txt  | awk '{print $5}')
		echo "${BASE1},${LIBNAME},${outputReads},with_primers,${sample}" >>  "${OUTPUT_SUMMARY}"

	done	
	
	
### End of the other pipeline





	 done


else #In case you already demultiplexed your samples, then cp the files you need
	cp "${DEMULT_OUTPUT}"/sample_trans.tmp "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/barcodes.fasta "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/summary.csv "${OUTPUT_DIR}"
	cp "${DEMULT_OUTPUT}"/pcr_primers.fasta "${OUTPUT_DIR}"

	DEMULT_DIR="${DEMULT_OUTPUT}"/demultiplexed

fi
	

 module unload cutadapt/4.1
 module unload miniconda/3.11



if [[ "${SEARCH_ASVs}" = "YES" ]]; then

	module load R/4.3.1

	Rscript "${SCRIPT_DIR}"/r/dada2_cluster.r "${OUTPUT_DIR}" "${DEMULT_DIR}" "${USE_HASH}"
fi
