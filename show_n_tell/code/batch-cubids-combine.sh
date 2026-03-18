#!/bin/bash

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate data-org

WRKDIR=$1
BASE=$2

echo python -u $WRKDIR/run-clin-dataorg.py $WRKDIR combine -b $BASE
python -u $WRKDIR/run-clin-dataorg.py $WRKDIR combine -b $BASE