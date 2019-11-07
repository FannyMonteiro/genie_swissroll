#!/bin/bash
#
# Script to restart ensemble run that is half-way through
# by Greg Colbourn (g.colbourn@uea.ac.uk)
#
# [1] set number of parts (e.g. 3 for 2-stage spin-up then experiment)
if [ -z "$1" ]
  then
    echo "1st argument is number of parts to experiment (e.g. 3 for 2-stage spin-up then experiment)"
    exit 60
  else
    NPARTS="$1"
fi
# [2] set number of years in each chunk
if
 [ -z "$2" ]; then
    echo "Usage: '$2' 2nd parameter must be the number of years in each chunk of model run"
    exit 65
  else
    MAXYEARS="$2"
fi
# [3] set the crash tolerance - min number of seconds allowed between resumits before job is killed
if
 [ -z "$3" ]; then
    echo "Usage: '$3' 8th parameter must be the crash tolerance - min number of seconds allowed between resumits before job is killed"
    exit 65
  else
    MINJOBTIME="$3"
fi
#
date
#
#ROOT=$PWD/../../..
#for purposes of keeping output pathnames short for netCDF:
ROOT=`echo "$PWD/../../.." | sed 's|/cgenie/genie-main/scripts/../../..|/cgenie/..|'`
SCRIPTSDIR=$PWD
EXECUTEROOT=$ROOT/cgenie_archive/to_resub
ARCHIVE=$ROOT/cgenie_archive
OUTPUTDIR=cgenie_output
PREFIX=cgenie
NEWPART=0
#
if [ -e resub_cgenie_myr_multipart_exe.sh ]
   then
   rm -f resub_cgenie_myr_multipart_exe.sh
fi
cp -f $SCRIPTSDIR/qsub_base.sh $SCRIPTSDIR/resub_cgenie_myr_multipart_exe.sh
echo "pwd" >> resub_cgenie_myr_multipart_exe.sh
echo ". /etc/profile.d/modules.sh" >> resub_cgenie_myr_multipart_exe.sh
echo "module add shared sge" >> resub_cgenie_myr_multipart_exe.sh
#
if [ -e fullnames ]	
	then   		
	rm -f fullnames
fi
if [ -e modelids ]	
	then   		
	rm -f modelids
fi
if [ -e baselines ]	
	then   		
	rm -f baselines
fi
if [ -e runids ]	
	then   		
	rm -f runids
fi	
cd $EXECUTEROOT
ls -l | grep $PREFIX | awk '	{print ""$9""}' > $SCRIPTSDIR/fullnames
sed 's\[.]\ \1' $SCRIPTSDIR/fullnames | sed 's\[.]\ \1' | awk '{print ""$1""}' > $SCRIPTSDIR/modelids
sed 's\[.]\ \1' $SCRIPTSDIR/fullnames | sed 's\[.]\ \1' | awk '{print ""$2""}' > $SCRIPTSDIR/baselines
sed 's\[.]\ \1' $SCRIPTSDIR/fullnames | sed 's\[.]\ \1' | awk '{print ""$3""}' > $SCRIPTSDIR/runids
J=0
NRUNS=`cat -n $SCRIPTSDIR/runids | tail -n1 | gawk '{print $1}'`
#NRUNS=`expr $NRUNS + 1`
while [ $J -lt $NRUNS ]
do
	J=`expr $J + 1`
	MODELID=`head -$J $SCRIPTSDIR/modelids | tail -1`
	echo "modelid = $MODELID"
	BASELINE=`head -$J $SCRIPTSDIR/baselines | tail -1`
	echo "baseline = $BASELINE"
    RUNID=`head -$J $SCRIPTSDIR/runids | tail -1`
	echo "runid = $RUNID"
	
	FULLNAME=`head -$J $SCRIPTSDIR/fullnames | tail -1`
    echo "$EXECUTEROOT/$FULLNAME"
	cd $EXECUTEROOT/$FULLNAME
	
	# clear pre-existing files	
	if [ -e tarballs ]
		then
   		rm -f tarballs
	fi
	
	# Only tarballs bigger than a set size [1000000] are used - smaller ones signify a model crash
	ls -l --sort=time |
	awk '	{if (NR > 1)
		{if ($5 > 1000000)
		{print ""$9""}}}' > $SCRIPTSDIR/tarballs
	FILENAME=`head -1 $SCRIPTSDIR/tarballs`
	echo "$FILENAME" > $SCRIPTSDIR/filename
	echo "FILENAME = $FILENAME"
	cd $SCRIPTSDIR
    #follow shortening done in runcgenie.sh
    source nameshortening.sh $MODELID.$RUNID 1
    NEWID=$SHORTNAME
    echo "NEWID = $NEWID"
	K=`sed "s/$NEWID.//g
	    s/.tar.gz//g
	     s\[.]\ \1 " $SCRIPTSDIR/filename | awk '{print ""$1""}'`
	echo "K = $K"
	#
	cd $EXECUTEROOT/$FULLNAME
	chmod 755 *
	if [ ! -e "$FILENAME" ]    # Check if file exists.
	then
		continue           # On to next.
	fi
	#unzip tarball
	gunzip -c $FILENAME | tar xvf - 
	LASTOUTPUT=`tail -1 $NEWID"_"$K/biogem/biogem_series_atm_pCO2.dat | awk '{print ""$1""}' | sed "s/[.]...//"`
	echo "last output = $LASTOUTPUT"
	NCHUNKS=`expr $LASTOUTPUT / $MAXYEARS`
    echo "number of chunks = $NCHUNKS"
    STARTYEAR=`expr $MAXYEARS \* \( $NCHUNKS + 1 \)`
    echo "start year = $STARTYEAR"
	#move restart files to cgenie_output
	echo "results folder = $NEWID"_"$K"
          cd $ROOT/$OUTPUTDIR
          chmod 755 *
	rm -rf $NEWID"_"$K
	cd $EXECUTEROOT/$FULLNAME
	mv -f $NEWID"_"$K $ROOT/$OUTPUTDIR/$NEWID"_"$K
	# write to the generated script
	cd $SCRIPTSDIR
	echo "/bin/bash qsub_cgenie_myr_multipart.sh $MODELID $BASELINE $RUNID $NPARTS $K $NEWPART $MAXYEARS $STARTYEAR $MINJOBTIME" >> resub_cgenie_myr_multipart_exe.sh
          echo "sleep 3" >> resub_cgenie_myr_multipart_exe.sh
done	
# run the generated script
cd $SCRIPTSDIR
chmod 755 resub_cgenie_myr_multipart_exe.sh
qsub resub_cgenie_myr_multipart_exe.sh
echo "all `expr $NRUNS` resubmitted"
# move archive folders back to right place
mv $EXECUTEROOT/* $EXECUTEROOT/../
