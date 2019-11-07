#for purposes of keeping output pathnames short for Mathematica
if [ -e shortname ]
   then
   rm -f shortname
fi
pwd
# for output directory names
if [ $2 -eq 1 ]
   then
   SHORTNAME=`echo "$1" | sed 's/genie_eb_go_gs_ac_bg_sg_//'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/ensemble_//'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/fullCC_//'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/preindustrial_fullCC_//'`
fi
# for script names
if [ $2 -eq 2 ]
   then
   SHORTNAME=`echo "$1" | sed 's/genie_eb_go_gs_ac_bg_sg_rg.//'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/el.//'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/ensemble_/e/'`
   SHORTNAME=`echo "$SHORTNAME" | sed 's/preindustrial_fullCC_//'`
fi
echo "name shortened: $1 -> $SHORTNAME"
