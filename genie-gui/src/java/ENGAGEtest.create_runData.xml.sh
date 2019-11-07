#! /bin/sh

source ../../../genie-main/user.sh
OUTFILE=../../test/runData.xml
echo "<run>" > $OUTFILE
echo "	<name>ENGAGEtest</name>" >> $OUTFILE
echo "	<outputDir>${OUTROOT}</outputDir>" >> $OUTFILE
echo "	<runDir>ENGAGEtest</runDir>" >> $OUTFILE
echo "	<initialDir>ENGAGEtest</initialDir>" >> $OUTFILE
echo "	<program>./genie_wrapper.sh</program>" >> $OUTFILE
echo "	<definitionXMLFile>definition.xml</definitionXMLFile>" >> $OUTFILE
echo "	<jobXMLFile>ENGAGEtest.xml</jobXMLFile>" >> $OUTFILE
echo "	<graphVarFile>graphVarData.xml</graphVarFile>" >> $OUTFILE
echo "	<mapVarFile>mapVarData.xml</mapVarFile>" >> $OUTFILE
echo "	<RUNTIME_ROOT>${RUNTIME_ROOT}</RUNTIME_ROOT>" >> $OUTFILE
echo "	<RUNTIME_OUTDIR>${RUNTIME_OUTDIR}</RUNTIME_OUTDIR>" >> $OUTFILE
echo "</run>" >> $OUTFILE
echo "<run>" >> $OUTFILE
echo "	<name>ENGAGEtest (Windows version)</name>" >> $OUTFILE
echo "	<outputDir>${OUTROOT}</outputDir>" >> $OUTFILE
echo "	<runDir>ENGAGEtest</runDir>" >> $OUTFILE
echo "	<initialDir>ENGAGEtest</initialDir>" >> $OUTFILE
echo "	<program>genie</program>" >> $OUTFILE
echo "	<definitionXMLFile>definition.xml</definitionXMLFile>" >> $OUTFILE
echo "	<jobXMLFile>ENGAGEtest.xml</jobXMLFile>" >> $OUTFILE
echo "	<graphVarFile>graphVarData.xml</graphVarFile>" >> $OUTFILE
echo "	<mapVarFile>mapVarData.xml</mapVarFile>" >> $OUTFILE
echo "	<RUNTIME_ROOT>${RUNTIME_ROOT}</RUNTIME_ROOT>" >> $OUTFILE
echo "	<RUNTIME_OUTDIR>${RUNTIME_OUTDIR}</RUNTIME_OUTDIR>" >> $OUTFILE
echo "</run>" >> $OUTFILE
cp ../../test/genie_wrapper.sh ${OUTROOT}/ENGAGEtest/
chmod u+x ${OUTROOT}/ENGAGEtest/genie_wrapper.sh
