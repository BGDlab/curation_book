#!/bin/bash

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate data-org

WRKDIR=$1
NAME=$2


echo "============ PREP ============"
echo python -u $WRKDIR/run-clin-dataorg.py $WRKDIR cubids-prep -n $NAME
python -u $WRKDIR/run-clin-dataorg.py $WRKDIR cubids-prep -n $NAME

if [ $? -eq 0 ] ; then
  echo "============ APPLY ============"
  echo python -u $WRKDIR/run-clin-dataorg.py $WRKDIR cubids-apply -n $NAME
  python -u $WRKDIR/run-clin-dataorg.py $WRKDIR cubids-apply -n $NAME

  if [ $? -eq 0 ] ; then
    echo "-------------"
  else
     echo FAILED TO LAUNCH CUBIDS APPLY; exit 1
  fi

else
  echo FAILED TO LAUNCH CUBIDS PREP; exit 1
fi