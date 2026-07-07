#!/bin/bash

#Script Loic VUILLERMET (UR ESCAPE)
#Creation data: 24 april 2025
#Last modification: 24 april 2025

#################################################################################################################
#                                                                                                               #
#       Needed: *_sorted.bam                                                                                    #
#                                                                                                               #
#       Remember: load the module/download the Apptainer Image: qualimap or having it at the good path          #
#                                                                                                               #
#################################################################################################################

# Add Qualimap to PATH (modifie ce chemin selon ton installation locale)
export PATH="/home/adm-loc/Documents/apps/qualimap:$PATH"

THREADS=12

# Verify the qualimap is available
if ! command -v qualimap &> /dev/null
then
    echo "Error: Qualimap not found in PATH!"
    exit 1
fi

# Prepare the output file
OUTPUT_TSV="qualimap_summary.tsv"
echo -e "sample\ttotal_read\tmapped_read\tpercentage_mapped\tmean_coverage\t1X\t5X\t10X" > "${OUTPUT_TSV}"

# Make a loop on all the *_sorted.bam
for BAM_FILE in *.sorted.bam
do
    ACCNUM=$(basename "${BAM_FILE}" .sorted.bam)
    OUTPUT_DIR="${ACCNUM}_qualimap"
    RESULT_FILE="${OUTPUT_DIR}/genome_results.txt"
    COVERAGE_FILE="${OUTPUT_DIR}/raw_data_qualimapReport/genome_fraction_coverage.txt"
    mkdir -p "${OUTPUT_DIR}"

    # Run the qualimap if necessary
    if [[ ! -s "${RESULT_FILE}" || "${BAM_FILE}" -nt "${RESULT_FILE}" ]]
    then
        echo "$(date '+%d-%m-%Y %H:%M:%S')- Running Qualimap for ${ACCNUM}..." >> log.txt
        qualimap bamqc \
        -bam "${BAM_FILE}" \
        -outdir "${OUTPUT_DIR}" \
        -outformat HTML \
        -nt ${THREADS}
    else
        echo "$(date '+%d-%m-%Y %H:%M:%S')- SKIP: ${ACCNUM} (up-to-date)." >> log.txt
    fi

    # Extract data (even if qualimap was not rerun)
    if [[ -s "${RESULT_FILE}" ]]
    then
        # Extract data from genome_results.txt
        TOTAL_READ=$(grep "number of reads" "${RESULT_FILE}" | awk -F'=' '{print $2}' | tr -d ' ,')
        MAPPED_READ=$(grep "number of mapped reads" "${RESULT_FILE}" | awk -F'=' '{print $2}' | awk '{print $1}' | tr -d ' ,')
        PERCENT_MAPPED=$(grep "number of mapped reads" "${RESULT_FILE}" | awk -F'(' '{print $2}' | awk -F'%' '{print $1}' | tr -d ' ')
        MEAN_COV=$(grep "mean coverage" "${RESULT_FILE}" | awk -F'=' '{print $2}' | tr -d ' ')

        # Extract data from genome_fraction_coverage.txt
        if [[ -f "${COVERAGE_FILE}" ]]
        then
            COV_1X=$(grep -w "1.0" "${COVERAGE_FILE}" | awk '{print $2}')
            COV_5X=$(grep -w "5.0" "${COVERAGE_FILE}" | awk '{print $2}')
            COV_10X=$(grep -w "10.0" "${COVERAGE_FILE}" | awk '{print $2}')
        else
            COV_1X="NA"
            COV_5X="NA"
            COV_10X="NA"
            echo "$(date '+%d-%m-%Y %H:%M:%S')- Error: No '${COVERAGE_FILE}'." >> log.txt
        fi

        echo -e "${ACCNUM}\t${TOTAL_READ}\t${MAPPED_READ}\t${PERCENT_MAPPED}\t${MEAN_COV}\t${COV_1X}\t${COV_5X}\t${COV_10X}" >> "${OUTPUT_TSV}"
        echo "$(date '+%d-%m-%Y %H:%M:%S')- Summary added for ${ACCNUM}..." >> log.txt
    else
        echo "$(date '+%d-%m-%Y %H:%M:%S')- WARNING: Result file missing or empty for ${ACCNUM}!" >> log.txt
    fi
done

# Sort the file (alphabetical) by sample
HEADER=$(head -n 1 "${OUTPUT_TSV}")
tail -n +2 "${OUTPUT_TSV}" | sort -k1,1 > "${OUTPUT_TSV}.tmp"
echo "${HEADER}" > "${OUTPUT_TSV}"
cat "${OUTPUT_TSV}.tmp" >> "${OUTPUT_TSV}"
rm "${OUTPUT_TSV}.tmp"
echo "$(date '+%d-%m-%Y %H:%M:%S')- JOB DONE." >> log.txt
echo ""
