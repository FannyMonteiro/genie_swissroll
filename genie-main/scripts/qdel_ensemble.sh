#!/bin/bash
#
# Script to quit ensemble jobs on cluster
# by Greg Colbourn (g.colbourn@uea.ac.uk)
#
date
#
# argument is ensemble number (e.g. e02a)
#
SCRIPTSDIR=$PWD
qstat | grep $1 |
awk '	{print ""$1""}' > $SCRIPTSDIR/jobids
NRUNS=`cat -n $SCRIPTSDIR/jobids | tail -n1 | gawk '{print $1}'`
echo "nruns = $NRUNS"
J=0
while [ $J -lt $NRUNS ]
do
	J=`expr $J + 1`
	RUNID=`head -$J $SCRIPTSDIR/jobids | tail -1`
         echo "$RUNID"
         qdel $RUNID
done
