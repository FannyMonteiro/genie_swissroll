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
#module add shared sge
#module add pgi/7.0.7
#module add netcdf-4.0
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
#ROOT=$PWD/../../..
#for purposes of keeping output pathnames short for netCDF:
ROOT=`echo "$PWD/../../.." | sed 's|/cgenie/genie-main/scripts/../../..|/cgenie/..|'`
PATCHES=$ROOT/cgenie/genie-userconfigs/_tmp.greg
SCRIPTSDIR=$PWD
ARCDIR=$ROOT/cgenie_archive
RUNLOG=$ROOT/cgenie_log
#
# add baseline config to runid (make sure it's done before so that anything specified in runid config overwrites baseline)
cat $PATCHES/$BASELINE"_0" > $PATCHES/$RUNID"_0"
cat $PATCHES/$RUNID >> $PATCHES/$RUNID"_0"
# make if loop for experinemt parts (e.g. to go from 1=spin-up1 to 2=spin-up2 to 3=experiment)
#
if [ $K -le $NPARTS ]
  then
  #Keep old MAXYEARS in case altered later on.
  MAXYEARSOLD=$MAXYEARS
  if [ $K -le 2 ]
     then
         STARTYEAR=0
     else
         STARTYEAR=$J
  fi
  # output directory for current part of experiment
  #OUTDIR=$ROOT/cgenie_output/$MODELID.$BASELINE.$K.$RUNID
  #for purposes of keeping output pathnames short for Mathematica:
  source nameshortening.sh $MODELID.$RUNID"_"$K 1
  NEWOUTDIR=$SHORTNAME
  echo "short outdir name = $NEWOUTDIR"
  OUTDIR=$ROOT/cgenie_output/$NEWOUTDIR
  #
  #output directory for previous part of experiment
  # now done so only a maximum 2 part spin-up
  O=`expr $K - 1`
  if [ $K -le 2 ]
     then
         L=`expr $K - 1`
     else
         L=2
  fi
  #OUTDIRPREV=$ROOT/cgenie_output/$MODELID.$BASELINE.$L.$RUNID
  #for purposes of keeping output pathnames short for Mathematica:
  source nameshortening.sh $MODELID.$RUNID"_"$L 1
  NEWOUTDIRPREV=$SHORTNAME
  echo "short old outdir name = $NEWOUTDIRPREV"
  OUTDIRPREV=$ROOT/cgenie_output/$NEWOUTDIRPREV
  #OUTDIRPREV is different for purposes of spin-up and restart for runs with history/future (i.e. 4 or more parts)
  #have OUTDIRPREV for spin-up purposes (as before) and OUTDIRPREVR for restart purposes
  source nameshortening.sh $MODELID.$RUNID"_"$O 1
  NEWOUTDIRPREV=$SHORTNAME
  echo "short old outdir name = $NEWOUTDIRPREV"
  OUTDIRPREVR=$ROOT/cgenie_output/$NEWOUTDIRPREV  
  #
  # add part-specific config to runid config
  cat $PATCHES/$RUNID"_0" > $PATCHES/$RUNID"_"$K
  cat $PATCHES/$BASELINE"_"$K >> $PATCHES/$RUNID"_"$K
  #
  # perform tasks specific to experiment part
  source $PATCHES/$BASELINE"_"$K".sh"
  sleep 3
  echo "run length: $RUNLENGTH"
  #
  echo ""    
  echo "config:"
  cat $PATCHES/$RUNID"_"$K
  echo "" 
  echo "year (J): $J"
  echo "" 
  echo "NEWPART: $NEWPART"
#
# make if loop for repeat running until runlength is reached
#
if [ $J -lt $RUNLENGTH ]
  then
	# time the gap between script submissions - if too short then model has crashed - see 'after' above
    before="$(date +%s)"
    # write in start year
	echo ">> Replacing simulation start year in goin with $J and setting ents_iwstp for runlength"
    # if RUNLENGTH - J < MAXYEARS then change runlength given to runcgenie.sh
    YEARSLEFT=`expr $RUNLENGTH - $J`
    echo "there are $YEARSLEFT years left to run for this part of the experiment"
    if [ $YEARSLEFT -lt $MAXYEARS ]
       then 
       MAXYEARS=$YEARSLEFT
       echo "Maxyears has been changed to $MAXYEARS"
    fi  
	IWSTP=`expr $MAXYEARS \* 100`
  	if [ -e $PATCHES/$RUNID"_"$K ]
          then
              echo "" >> $PATCHES/$RUNID"_"$K
	          echo "bg_par_misc_t_start=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "rg_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "sg_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "ents_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "ents_iwstp=$IWSTP" >> $PATCHES/$RUNID"_"$K
          else
              touch $PATCHES/$RUNID"_"$K
              echo "" >> $PATCHES/$RUNID"_"$K
              echo "ents_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "sg_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "rg_start_year=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "bg_par_misc_t_start=$J.0" >> $PATCHES/$RUNID"_"$K
              echo "ents_iwstp=$IWSTP" >> $PATCHES/$RUNID"_"$K
          fi
    #   
    # run the model
	# (run from restart unless first run)
    if [ "$NEWPART" = "1" ]
	then
		cd $SCRIPTSDIR
  		echo "running runcgenie script..."
  		if [ $K -gt 1 ]
  		# if continuing (parts 2 or 3 or 4 or 5) select restart id:
  		then
            echo "./runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS $OUTDIRPREVR"
  			/bin/bash runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS $OUTDIRPREVR
  		# otherwise run from scratch (part 1):
 		else
 		    echo "./runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS"
 			/bin/bash runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS
  		fi
        # create archive direcory
        mkdir $ARCDIR/$MODELID.$BASELINE.$RUNID
        NEWPART=0
	else	
		cd $SCRIPTSDIR
  		echo "running runcgenie script..."
            echo "./runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS $OUTDIR"
 			/bin/bash runcgenie.sh $MODELID $RUNID"_"$K $MAXYEARS $OUTDIR
	fi
	#
	# increment J by $MAXYEARS
	J=`expr $J + $MAXYEARS`
	# adjust newid for submission script name (make condensed so can see in job list)
    source nameshortening.sh $RUNID"_"$K 2
    NEWRUNID=$SHORTNAME
    echo "short runid = $NEWRUNID"
	#NEWJ=`echo "$J" | sed s/000//`
	#NEWID=$NEWRUNID"_"$NEWJ
	NEWID=$NEWRUNID"_"$J
    echo "NEWID = $NEWID"
	#
	# move output to archive directory	
	cd $ARCDIR/fresh
	pwd
	mv -v $NEWOUTDIR.tar.gz $ARCDIR/$MODELID.$BASELINE.$RUNID/$NEWOUTDIR.$J.tar.gz
	cd $SCRIPTSDIR
	pwd
	#
	# change J in $RUNLOG/$NEWID.sh script and resubmit job
#	cp -f qsub_base.sh $RUNLOG/$NEWID.sh
#	echo "cd $SCRIPTSDIR" >> $RUNLOG/$NEWID.sh
#	echo "pwd" >> $RUNLOG/$NEWID.sh
#	echo "/bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $J $MINJOBTIME" >> $RUNLOG/$NEWID.sh
#	cat $RUNLOG/$NEWID.sh
    after="$(date +%s)"
    # if elapsed seconds between script submissions is too short then model has crashed so exit
    elapsed_seconds="$(expr $after - $before)"
    if [ $elapsed_seconds -lt $MINJOBTIME ]
       then 
           echo "model has most likely crashed - only $elapsed_seconds seconds between successive script submits (min allowed is $MINJOBTIME)!"
           exit
       else
       	# delete old archive tarball here instead so if model crashes still have files to restart
       	OLDJ=`expr $J - $MAXYEARS`
	    rm -f $ARCDIR/$MODELID.$BASELINE.$RUNID/$NEWOUTDIR.$OLDJ.tar.gz
    fi
#     quotacheck > quotacheckfile
#     USED="$(cat quotacheckfile | grep envhome -1 | tail -1 | gawk '{printf ("%.0f", $3)}')"
#     QUOTA="$(cat quotacheckfile | grep envhome -1 | tail -1 | gawk '{printf ("%.0f", $5)}')"
#     LEFT="$(expr $QUOTA - $USED)"
#     echo "disk space = $USED/$QUOTA; $LEFT left"
#     if [ $LEFT -lt 10 ]
#     then
#         echo "remaing quota < 10 Gig: STOPPING"
#         exit
#     fi
#	qsub $RUNLOG/$NEWID.sh
     # change script to run here rather than as job because of nodes not being "submit hosts"
     cd $SCRIPTSDIR
     pwd
 	 echo "/bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $J $MINJOBTIME"
     /bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $J $MINJOBTIME

  else
  	 echo ""
	 echo "Part $K is finished..."
	 echo ""
	 echo ""
	 # increment K by 1
     K=`expr $K + 1`
     NEWPART=1
     echo "NEWPART: $NEWPART"
     # adjust newid for submission script name (make condensed so can see in job list)
     source nameshortening.sh $RUNID"_"$K 2
     NEWRUNID=$SHORTNAME
     echo "short runid = $NEWRUNID"
	 NEWID=$NEWRUNID"_$STARTYEAR"
     # change K in $RUNLOG/$NEWID.sh script and resubmit job
# 	 cp -f qsub_base.sh $RUNLOG/$NEWID.sh
# 	 echo "cd $SCRIPTSDIR" >> $RUNLOG/$NEWID.sh
# 	 echo "pwd" >> $RUNLOG/$NEWID.sh
# 	 echo "/bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $STARTYEAR $MINJOBTIME" >> $RUNLOG/$NEWID.sh
# 	 cat $RUNLOG/$NEWID.sh
#      qsub $RUNLOG/$NEWID.sh
     # change script to run here rather than as job because of nodes not being "submit hosts"
     cd $SCRIPTSDIR
     pwd
 	 echo "/bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $STARTYEAR $MINJOBTIME"
     /bin/bash cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARSOLD $STARTYEAR $MINJOBTIME
     exit
  fi
else
	echo "Cleaning up run logs..."
#	cp -f $SCRIPTSDIR/qsub_base.sh $SCRIPTSDIR/qsub_sort_runlog.sh
#    echo "pwd" >> $SCRIPTSDIR/qsub_sort_runlog.sh
#    echo "/bin/bash sort_runlog_date.sh $RUNLOG" >> $SCRIPTSDIR/qsub_sort_runlog.sh
#	cat $SCRIPTSDIR/qsub_sort_runlog.sh
#	qsub $SCRIPTSDIR/qsub_sort_runlog.sh
	pwd
	echo "/bin/bash sort_runlog_date.sh $RUNLOG"
	/bin/bash sort_runlog_date.sh $RUNLOG
    exit
fi

