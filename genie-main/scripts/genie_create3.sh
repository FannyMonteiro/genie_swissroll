#!/bin/bash
#
# Script to create a new instance of genie for running;
# the first argument supplied is the name of the directory it is contained in.
# the second argument supplied is the top level directory containing the first argument's directory
# by Greg Colbourn (g.colbourn@uea.ac.uk)
#
date
# make directories
cd $2
mkdir $1
cd $1
# checkout from svn repository
svn checkout http://source.ggy.bris.ac.uk/subversion/genie/branches/cgenie --username greg-colbourn cgenie
# make needed directories
mkdir cgenie_archive
mkdir cgenie_log
cd genie/genie-main
# edit user.* files to give correct paths
sed s/genie_dev/$1/g < user.mak > user2.mak
mv user2.mak user.mak
sed s/genie_dev/$1/g < user.sh > user2.sh
mv user2.sh user.sh
# make and run tests
make cleanall
make assumedgood
make test
make testebgogs
make testents
make testbiogem
# do test run
cd ../genie-tools/runscripts
./rungenie.sh genie_eb_go_gs_ac_bg_sg_rg_el rokgem_test 10 


