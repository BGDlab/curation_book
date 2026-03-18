#!/bin/bash

# This script launches an array of jobs to convert a sourcedata directory to a BIDS directory.


WRKDIR=$1
DATASET=$2
PC=$3


BIDS="${DATASET}/BIDS"
SRC="${DATASET}/sourcedata"
CONFIG=$4
SCRIPT="${WRKDIR}/data-bidsify.py"


line="`sed -n ${SLURM_ARRAY_TASK_ID}p ${CONFIG}`"

subj=$(echo $line | cut -d"," -f1)      # e.g. HMUXLIVN
sess=$(echo $line | cut -d',' -f4)      # e.g. 65443710698procId004206ageDays or 65443710698procId004206ageDays 168153109945procId003840ageDays
dirs=$(echo $line | cut -d',' -f5)      # e.g. $DICOM_SRC/HMUXLIVN_65443710698_4206 or $DICOM_SRC/HMUXLIVN_65443710698_4206 $DICOM_SRC/HMUXLIVN_168153109945_3840


######## RUN BIDSIFICATION WITH HEUDICONV ########
if [ $PC == 'False' ]; then
  
  #DCMS=$5

  echo python -u $SCRIPT $WRKDIR -d $dirs -i $SRC -o $BIDS -s $subj -ss $sess
  python -u $SCRIPT $WRKDIR -d $dirs -i $SRC -o $BIDS -s $subj -ss $sess
    



######## DO POST CUBIDS UPDATES ########
else

  MAP=$5

  echo python -u $SCRIPT $WRKDIR -o $BIDS -s $subj -ss $sess -pc -mm $MAP
  python -u $SCRIPT $WRKDIR -o $BIDS -s $subj -ss $sess -pc -mm $MAP

fi