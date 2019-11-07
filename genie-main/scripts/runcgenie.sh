#!/bin/bash
#
#####################################################################
### SCIPT TO META-CONFIGURE AND RUN GENIE: CLIMATE + CARBON cycle ###
#####################################################################
#
# (0) USER OPTIONS
# ----------------
#
HOME=/home/lobster/cgenie
#
# (1) GET PASSED PARAMETERS
# -------------------------
# [1] base configuration ID
if [ -z "$1" ]; then
    echo "Usage: '$1' 1st parameter must be the config ID"
    exit 65
  else
    MODELID="$1"
fi
# [2] set run ID (input run ID (= configuration patch file name))
if [ -z "$2" ]; then
    echo "Usage: '$2' 2nd parameter must be the run ID"
    exit 65
  else
    RUNID="$2"
fi
# [3] set run duration
if [ -z "$3" ]; then
    echo "Usage: '$3' 3rd parameter must be the run length (years)"
    exit 65
  else
    RUNLENGTH="$3"
fi
# [4] restart path (optional)
if [ -n "$4" ]; then
  RESTARTPATH="$4"
fi
#
# (2) SET LOCAL FILE AND DIRECTORY NAMES
# --------------------------------------
# ensure rocks can find xsltproc
export PATH=$PATH:/opt/rocks/bin:/usr/bin
export PATH=$PATH:/share/apps/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/lib
#
#ROOT=$PWD/../../..
#for purposes of keeping output pathnames short for netCDF:
ROOT=`echo "$PWD/../../.." | sed 's|/cgenie/genie-main/scripts/../../..|/cgenie/..|'`
GOINDIR=$ROOT/cgenie/genie-userconfigs/_tmp.greg
GOIN=$GOINDIR/$RUNID
#OUTPUTPATH=genie_output/$MODELID.$RUNID
#for purposes of keeping output pathnames short for Mathematica:
source nameshortening.sh $MODELID.$RUNID 1
NEWID=$SHORTNAME
OUTPUTPATH=$ROOT/cgenie_output/$NEWID
echo "output path: $OUTPUTPATH"
#$OUTPUTPATH=out/$NEWID
CONFIGPATH=$ROOT/cgenie/genie-main/configs
CONFIGNAME=$MODELID"."$RUNID".config"
BINARYPATH=$ROOT/cgenie/genie-main
RESTARTNAME="rst.1"
if test -e $CONFIGPATH/$MODELID".config"
then
    echo "   Experiment configuration: "
    echo $CONFIGPATH/$MODELID".config"
    echo " found :)"
else
    echo "   Experiment configuration: "
    echo $CONFIGPATH/$MODELID".config"
    echo " CANNOT be found :("
    exit 1
fi
if test -d $GOINDIR
then
    echo "   Experiment configuration (namelist changes) directory: "
    echo $GOINDIR
    echo " found :)"
else
    echo "   Experiment configuration (namelist changes) directory: "
    echo $GOINDIR
    echo " CANNOT be found :("
    exit 1
fi
if test -e $GOINDIR/$RUNID
then
    echo "   Experiment configuration patch file: "
    echo $GOINDIR/$RUNID
    echo " found :)"
else
    echo "   Experiment configuration patch file: "
    echo $GOINDIR/$RUNID
    echo " CANNOT be found :("
    exit 1
fi
#
# (3) CREATE RUN CONFIG FILE
# --------------------------
# Copy template config file
cp -f $CONFIGPATH/$MODELID".config" $CONFIGPATH/$CONFIGNAME
# Set the experiment run name
#echo EXPID=$MODELID.$RUNID >> $CONFIGPATH/$CONFIGNAME
#for purposes of keeping output pathnames short for Mathematica:
echo EXPID=$NEWID >> $CONFIGPATH/$CONFIGNAME
#
# (4) SET MODEL TIME-STEPPING
# ---------------------------
# ************************************
# DEFINE NUMBER OF TIME-STEPS PER YEAR
let N_TIMESTEPS=96
# ************************************
# set BIOGEM run length
echo bg_par_misc_t_runtime=$RUNLENGTH >> $CONFIGPATH/$CONFIGNAME
# set SEDGEM sediment age
echo sg_par_misc_t_runtime=$RUNLENGTH >> $CONFIGPATH/$CONFIGNAME
# set overall GENIE run length
let stp=$RUNLENGTH*5*$N_TIMESTEPS
echo ma_koverall_total=$stp >> $CONFIGPATH/$CONFIGNAME
echo ma_dt_write=$stp >> $CONFIGPATH/$CONFIGNAME
# set climate model component restart frequency
let stp=$RUNLENGTH*$N_TIMESTEPS
echo ea_4=$stp >> $CONFIGPATH/$CONFIGNAME
echo go_4=$stp >> $CONFIGPATH/$CONFIGNAME
echo gs_4=$stp >> $CONFIGPATH/$CONFIGNAME
echo ents_iwstp=$stp >> $CONFIGPATH/$CONFIGNAME
# set 'health check' frequency [NOTE: a '+1' in effect disables this feature]
let stp=$RUNLENGTH*$N_TIMESTEPS
echo ea_3=$stp >> $CONFIGPATH/$CONFIGNAME
echo go_3=$stp >> $CONFIGPATH/$CONFIGNAME
echo gs_3=$stp >> $CONFIGPATH/$CONFIGNAME
echo ents_npstp=$stp >> $CONFIGPATH/$CONFIGNAME
# set 'time series' frequency [NOTE: a '+1' in effect disables this feature]
let stp=$RUNLENGTH*$N_TIMESTEPS+1
echo ea_5=$stp >> $CONFIGPATH/$CONFIGNAME
echo go_5=$stp >> $CONFIGPATH/$CONFIGNAME
echo gs_5=$stp >> $CONFIGPATH/$CONFIGNAME
echo ents_npstp=$stp >> $CONFIGPATH/$CONFIGNAME
# set 'an average' frequency [NOTE: a '+1' in effect disables this feature]
let stp=$RUNLENGTH*$N_TIMESTEPS+1
echo ea_6=$stp >> $CONFIGPATH/$CONFIGNAME
echo go_6=$stp >> $CONFIGPATH/$CONFIGNAME
echo gs_6=$stp >> $CONFIGPATH/$CONFIGNAME
echo ents_ianav=$stp >> $CONFIGPATH/$CONFIGNAME
# SET CLIMATE COMPONENTS TIME-STEPS PER YEAR
echo ea_9=$N_TIMESTEPS >> $CONFIGPATH/$CONFIGNAME
echo go_9=$N_TIMESTEPS >> $CONFIGPATH/$CONFIGNAME
echo gs_9=$N_TIMESTEPS >> $CONFIGPATH/$CONFIGNAME
#
# (5) SET CLIMATE MODEL RE-START FILE DETAILS
# -------------------------------------------
# Set default flags
# Set netCDF restart saving flag
echo ea_31=n >> $CONFIGPATH/$CONFIGNAME
echo go_19=n >> $CONFIGPATH/$CONFIGNAME
echo gs_14=n >> $CONFIGPATH/$CONFIGNAME
#echo ents_netout=n >> $CONFIGPATH/$CONFIGNAME
# Set ASCII restart output flag
echo ea_32=y >> $CONFIGPATH/$CONFIGNAME
echo go_20=y >> $CONFIGPATH/$CONFIGNAME
echo gs_15=y >> $CONFIGPATH/$CONFIGNAME
# Set ASCII restart number (i.e., output file string)
echo ea_29=rst >> $CONFIGPATH/$CONFIGNAME
echo go_17=rst >> $CONFIGPATH/$CONFIGNAME
echo gs_12=rst >> $CONFIGPATH/$CONFIGNAME
echo ents_out_name="rst" >> $CONFIGPATH/$CONFIGNAME
echo ents_restart_file="rst.sland" >> $CONFIGPATH/$CONFIGNAME
#
# (6) CONFIGURE USE OF RESTART
# -----------------------------
# Set continuing/new run flags
# => set restart input flags
# => disable netCDF restart input flag
# => set restart input number
# => copy restart files to data directory
if [ -n "$4" ]; then
  echo ">> Checking whether restart directory $RESTARTPATH exists ..."
  if test -d $RESTARTPATH
  then
    echo "   OK :)"
  else
      echo "   Restart directory $RESTARTPATH cannot be found"
    exit 1
  fi
  echo ea_7=c >> $CONFIGPATH/$CONFIGNAME
  echo go_7=c >> $CONFIGPATH/$CONFIGNAME
  echo gs_7=c >> $CONFIGPATH/$CONFIGNAME
  echo ents_restart=c >> $CONFIGPATH/$CONFIGNAME
  echo ac_ctrl_continuing=t >> $CONFIGPATH/$CONFIGNAME
  echo bg_ctrl_continuing=t >> $CONFIGPATH/$CONFIGNAME
  echo sg_ctrl_continuing=t >> $CONFIGPATH/$CONFIGNAME
  echo rg_ctrl_continuing=t >> $CONFIGPATH/$CONFIGNAME
  echo ea_30=n >> $CONFIGPATH/$CONFIGNAME
  echo go_18=n >> $CONFIGPATH/$CONFIGNAME
  echo gs_13=n >> $CONFIGPATH/$CONFIGNAME
  echo ents_netin=n >> $CONFIGPATH/$CONFIGNAME
  echo ea_35=$RESTARTNAME >> $CONFIGPATH/$CONFIGNAME
  echo go_23=$RESTARTNAME >> $CONFIGPATH/$CONFIGNAME
  echo gs_18=$RESTARTNAME >> $CONFIGPATH/$CONFIGNAME
  #echo ents_restart_file="spn.sland" >> $CONFIGPATH/$CONFIGNAME
  #cp $RESTARTPATH/ents/spn.sland cgenie/genie-ents/data/
  echo ea_rstdir_name=$RESTARTPATH"/embm" >> $CONFIGPATH/$CONFIGNAME
  echo go_rstdir_name=$RESTARTPATH"/goldstein" >> $CONFIGPATH/$CONFIGNAME
  echo gs_rstdir_name=$RESTARTPATH"/seaice" >> $CONFIGPATH/$CONFIGNAME
  echo ents_outdir_name=$OUTPUTPATH"/ents" >> $CONFIGPATH/$CONFIGNAME
  echo ents_dirnetout=$RESTARTPATH"/ents" >> $CONFIGPATH/$CONFIGNAME
  echo ents_rstdir_name=$RESTARTPATH"/ents" >> $CONFIGPATH/$CONFIGNAME
  echo gs_rstdir_name=$RESTARTPATH"/goldsteinseaice" >> $CONFIGPATH/$CONFIGNAME
  echo ac_par_rstdir_name=$RESTARTPATH"/atchem" >> $CONFIGPATH/$CONFIGNAME
  echo bg_par_rstdir_name=$RESTARTPATH"/biogem" >> $CONFIGPATH/$CONFIGNAME
  echo bg_par_outdir_name=$OUTPUTPATH"/biogem" >> $CONFIGPATH/$CONFIGNAME
  echo sg_par_rstdir_name=$RESTARTPATH"/sedgem" >> $CONFIGPATH/$CONFIGNAME
  echo sg_par_outdir_name=$OUTPUTPATH"/sedgem" >> $CONFIGPATH/$CONFIGNAME
  echo rg_par_rstdir_name=$RESTARTPATH"/rokgem" >> $CONFIGPATH/$CONFIGNAME
  echo rg_par_outdir_name=$OUTPUTPATH"/rokgem" >> $CONFIGPATH/$CONFIGNAME
else
  echo ea_7=n >> $CONFIGPATH/$CONFIGNAME
  echo go_7=n >> $CONFIGPATH/$CONFIGNAME
  echo gs_7=n >> $CONFIGPATH/$CONFIGNAME
  echo ents_restart=n >> $CONFIGPATH/$CONFIGNAME
  echo ac_ctrl_continuing=f >> $CONFIGPATH/$CONFIGNAME
  echo bg_ctrl_continuing=f >> $CONFIGPATH/$CONFIGNAME
  echo sg_ctrl_continuing=f >> $CONFIGPATH/$CONFIGNAME
  echo rg_ctrl_continuing=f >> $CONFIGPATH/$CONFIGNAME
fi
#
# (7) INCORPORATE RUN CONFIGURATION
# ---------------------------------
dos2unix $GOIN
cat $GOIN >> $CONFIGPATH/$CONFIGNAME
#
# (8) GO!
# -------
# Run model ...
echo ">> Here we go ..."
cd $BINARYPATH
./old_genie_example.job -f configs/$CONFIGNAME
# Archive and clean up
cd ../..
echo ">> Archiving ..."
# move config file to out
mv $CONFIGPATH/$CONFIGNAME $OUTPUTPATH/$CONFIGNAME
# archive results
cd cgenie_output
tar cfz ../cgenie_archive/fresh/$NEWID.tar.gz $NEWID
echo ">> All done - now go and play outside"
#
