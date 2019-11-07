#! /bin/sh

# ##############################################  #
# build_experiment.sh                             #
# script to generate genie experiment directory   #
# which is ready to go for use with ALADDIN GUI or#
# to run as normal. All input files, configs, svn #
# info etc recorded for tracibility (some of this #
# is done in the call to genie_example.job)       #
# authors Simon Mueller s.a.mueller@open.ac.uk    #
# and Martin Johnson martin.johnson@uea.ac.uk     #
# *********************************************** #

. ../genie-main/user.sh

# DEFAULT VALUES FOR OPTIONAL COMMAND LINE ARGUMENTS
MAP_VAR_DATA_FILE="FALSE"
GRAPH_VAR_DATA_FILE="FALSE"
CLEANALL="YES"


if [ -z $1 ]; then echo "Usage: ./build_experiment.sh -c <CONFIGURATION> -f
    <FILELIST> -m <MAPVARDATA FILE> -g <GRAPHVARDATA FILE> -n: don't make
    cleanall before running genie_example.job (not reccomended unless you know
    what you're doing!"; exit 1; fi
while getopts "nc:f:m:g:" opt; do
    case $opt in
	n) CLEANALL="NO" ;;
	c) EXPERIMENT_CONFIG=$OPTARG ;;
	f) EXPERIMENT_FILELIST=$OPTARG ;;
	m) MAP_VAR_DATA_FILE=$OPTARG ;;
	g) GRAPH_VAR_DATA_FILE=$OPTARG ;;
	*) echo "Usage: ./build_experiment.sh -c <CONFIGURATION> -f <FILELIST> -m <MAPVARDATA FILE> -g <GRAPHVARDATA FILE> -n: don't make cleanall before running genie_example.job (not reccomended unless you know what you're doing!"; exit 1;;
    esac
done

EXPERIMENT_ROOT=`pwd`
EXPERIMENT_NAME=`xsltproc getEXPID.xsl ${EXPERIMENT_CONFIG}`

echo "creating experiment directory..."
if [ ! -d ${EXPERIMENT_NAME} ]; then mkdir ${EXPERIMENT_NAME}; fi
echo "creating input data directory..."
if [ ! -d ${EXPERIMENT_NAME}/input ]; then mkdir ${EXPERIMENT_NAME}/input; fi
echo "creating aladdin directory..."
if [ ! -d ${EXPERIMENT_NAME}/aladdin ]; then mkdir ${EXPERIMENT_NAME}/aladdin; fi
echo "copying *varData files"

if [ ${GRAPH_VAR_DATA_FILE} != "FALSE" ] ; then 
    cp ${GRAPH_VAR_DATA_FILE} ${EXPERIMENT_NAME}/aladdin/graphVarData.xml
else
    cp ${CODEDIR}/genie-gui/graphVarData.xml ${EXPERIMENT_NAME}/aladdin/
fi 

if [ ${MAP_VAR_DATA_FILE} != "FALSE" ] ; then 
    cp ${MAP_VAR_DATA_FILE} ${EXPERIMENT_NAME}/aladdin/mapVarData.xml
else
    cp ${CODEDIR}/genie-gui/mapVarData.xml ${EXPERIMENT_NAME}/aladdin/
fi 

echo "copying config files..."
# copy file from source tree into experiment tree, but keep the files that
# already exist in the experiment tree (option 'k' for tar)..."
cd ${CODEDIR}
tar cf /tmp/$$.tar `cat ${EXPERIMENT_ROOT}/${EXPERIMENT_FILELIST}`
echo "extract" 
cd ${EXPERIMENT_ROOT}/${EXPERIMENT_NAME}/input
tar -kxf /tmp/$$.tar

echo "copying genie_wrapper.sh"
cp ${CODEDIR}/genie-gui/genie_wrapper.sh ${EXPERIMENT_ROOT}/${EXPERIMENT_NAME}
echo ${RUNTIME_OUTDIR}

echo "building experiment ${EXPERIMENT_NAME}..."
cp ${EXPERIMENT_ROOT}/${EXPERIMENT_CONFIG} ${CODEDIR}/genie-main/configs
cd ${CODEDIR}/genie-main
if [ ${CLEANALL} = "YES" ] ; then
    make cleanall
fi
./genie_example.job -x -f ${EXPERIMENT_CONFIG}

