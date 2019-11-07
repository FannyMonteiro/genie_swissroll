#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -j y
#$ -o ../../../cgenie_log/
pwd
/bin/bash sort_runlog_date.sh /home/lobster/cgenie/cgenie/../cgenie_log
