#!/bin/bash
#
# Script to sort directory (with many files in it)
# into subdirectories arranged by date (column 6 in ls -l).
# Can change 6 to other columns to sort by other things (e.g. 3 is owner).
# the path to the directory is given as an argument on the command line
# by Greg Colbourn - g.colbourn@uea.ac.uk
#
# note: #7 is for ESCLuster, should be #8 for cluster1
cd $1
if [ -e zz_sort_runlog_exe.sh ]
   then
   rm -f zz_sort_runlog_exe.sh
fi
ls -l  --time-style=+%F |
awk '	{if (NR > 1) 
	{print	"if [ -d "$6" ]",
		"\n\tthen",
		"\n\t\techo""",
		"\n\telse",
		"\n\t\tmkdir -p "$6,
		"\nfi",
		"\nif [ -f "$7" ]",
		"\n\tthen",
		"\nmv "$7" "$6"/",
		"\nfi"}}' > zz_sort_runlog_exe.sh
chmod 755 zz_sort_runlog_exe.sh
./zz_sort_runlog_exe.sh
