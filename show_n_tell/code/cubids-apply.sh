#!/usr/bin/bash

# Usage

# Run normally
#     sbatch cubids_apply.sh /path/to/working-dir 
#         /path/to/dataset-dir 
#         False

# Run batched
#     sbatch cubids_apply.sh /path/to/working-dir 
#         /path/to/dataset-dir 
#         /path/to/batch_config.csv


BGD=/mnt/isilon/bgdlab_processing
SIMG=/mnt/isilon/bgdlab_processing/code/containers/clin-dataorg-0.0.2.sif


WRKDIR=$1
DIR=$2
CONFIG=/app/scripts/config.yml
ARRAY=$3
BIDS="${DIR}/BIDS"

######## RUN NORMALLY ########
if [ $ARRAY == 'False' ]; then

    CUBIDS="${BIDS}/code/CuBIDS"
    TABLE="${CUBIDS}/run01_summary_edited.tsv"
    FILES="${CUBIDS}/run01_files.tsv"
    OUT='run02'


    echo running cubids apply
    echo singularity run -e -C -B $BGD $SIMG cubids apply --config $CONFIG $BIDS $TABLE $FILES $OUT
    singularity run -e -C -B $BGD $SIMG cubids apply --config $CONFIG $BIDS $TABLE $FILES $OUT

    # Look for the runXX_full_cmd.sh file and move it into the BIDS/code/CuBIDS directory to prevent 
    # it being overwritten by simulatenous runs on different datasets
    if [ -f "${WRKDIR}/_full_cmd.sh" ]; then
        echo find $WRKDIR -type f -name "_full_cmd.sh" -exec mv {} $CUBIDS \;
        find $WRKDIR -type f -name "_full_cmd.sh" -exec mv {} $CUBIDS \;
    fi


######## RUN BATCHED ########
else
    echo "config ${ARRAY}"
    line="`sed -n ${SLURM_ARRAY_TASK_ID}p ${ARRAY}`"

    CUBIDS="${DIR}/BIDS/code/CuBIDS/batch"

    TABLE=$(echo $line | cut -d"," -f2)
    FILES=$(echo $line | cut -d"," -f3)
    OUT=$(echo $line | cut -d"," -f4)


    echo running cubids apply
    echo singularity run -e -C -B $BGD $SIMG cubids apply --config $CONFIG $BIDS $TABLE $FILES $OUT
    singularity run -e -C -B $BGD $SIMG cubids apply --config $CONFIG $BIDS $TABLE $FILES $OUT

    # Look for the runXX_full_cmd.sh file and move it into the BIDS/code/CuBIDS directory to prevent 
    # it being overwritten by simulatenous runs on different datasets
    if [ -f "${WRKDIR}/_full_cmd.sh" ]; then
        echo find $WRKDIR -type f -name "_full_cmd.sh" -exec mv {} $CUBIDS \;
        find $WRKDIR -type f -name "_full_cmd.sh" -exec mv {} $CUBIDS \;
    fi

fi