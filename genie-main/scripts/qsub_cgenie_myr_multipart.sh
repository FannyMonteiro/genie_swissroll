#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -j y
#$ -o ../../../cgenie_log/
#
# Script to submit very long cgenie runs with sediments and weathering, including spin-ups - see GENIE HOW-TO document.
# by Greg Colbourn (g.colbourn@uea.ac.uk)
#
#. /etc/profile.d/modules.sh
# module add shared sge
#
date
#
# (1) GET PASSED PARAMETERS
# -------------------------
# [1] base configuration ID
if [ -z "$1" ]; then
    echo "Usage: '$1' 1st parameter must be the config ID e.g. cgenie_eb_go_gs_ac_bg"
    exit 65
  else
    MODELID="$1"
fi
# [2] set baseline config name
if [ -z "$2" ]; then
    echo "Usage: '$2' 2nd parameter must be the baseline config name, e.g. worbe_fulCC"
    exit 65
  else
    BASELINE="$2"
fi
# [3] set run ID (config patch file)
if [ -z "$3" ]; then
    echo "Usage: '$3' 3rd parameter must be the run ID e.g. worbe2_preindustrial_1"
    exit 65
  else
    RUNID="$3"
fi
# [4] set number of parts (e.g. 3 for 2-stage spin-up then experiment)
if [ -z "$4" ]
  then
    echo "4th argument is number of parts to experiment (e.g. 3 for 2-stage spin-up then experiment)"
    exit 60
  else
    NPARTS="$4"
fi
# [5] set experiment counter (to go from 1=spin-up1 to 2=spin-up2 to 3=experiment)
if [ -z "$5" ]
  then
    echo "5th argument is experiment index to icrement should be set to 1"
    exit 60
  else
    K="$5"
fi
# [6] Is it the start of a new part? (0=no, 1=yes - set to 1 initially)
if [ -z "$5" ]
  then
    echo "Is it the start of a new part? (0=no, 1=yes - set to 1 initially)"
    exit 60
  else
    NEWPART="$6"
fi
# [7] set number of years in each chunk
if
 [ -z "$7" ]; then
    echo "Usage: '$7' 6th parameter must be the number of years in each chunk of model run"
    exit 65
  else
    MAXYEARS="$7"
fi
# [8] set year counter
if [ -z "$8" ]
  then
    echo "8th argument is year index to icrement"
    exit 60
  else
    J="$8"
fi
# [9] set the crash tolerance - min number of seconds allowed between resumits before job is killed
if
 [ -z "$9" ]; then
    echo "Usage: '$9' 9th parameter must be the crash tolerance - min number of seconds allowed between resumits before job is killed"
    exit 65
  else
    MINJOBTIME="$9"
fi
#
echo "Arguments for cgenie_myr_multipart.sh: $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARS $J $MINJOBTIME"
#
# define root paths
SCRIPTSDIR=$PWD
#ROOT=$PWD/../../..
#for purposes of keeping output pathnames short for netCDF:
ROOT=`echo "$PWD/../../.." | sed 's|/cgenie/genie-main/scripts/../../..|/cgenie/..|'`
PATCHES=$ROOT/cgenie/genie-userconfigs/_tmp.greg
RUNLOG=$ROOT/cgenie_log
#
# submit cgenie_myr_multipart.sh script as qsub job
#
# adjust newid for submission script name (make condensed so can see in job list)
source nameshortening.sh $RUNID"_"$K 2
NEWRUNID=$SHORTNAME
echo "short runid = $NEWRUNID"
NEWID=$NEWRUNID"_$J"
# change K in $RUNLOG/$NEWID.sh script and resubmit job
cp -f qsub_base.sh $RUNLOG/$NEWID.sh
echo "cd $SCRIPTSDIR" >> $RUNLOG/$NEWID.sh
echo "pwd" >> $RUNLOG/$NEWID.sh
echo "/bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARS $J $MINJOBTIME" >> $RUNLOG/$NEWID.sh
echo " " >> $RUNLOG/$NEWID.sh
cat $RUNLOG/$NEWID.sh
qsub $RUNLOG/$NEWID.sh
#
echo "submitted $RUNLOG/$MODELID.$BASELINE.$K.$RUNID.sh"
#

