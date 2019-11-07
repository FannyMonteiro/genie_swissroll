quotacheck > quotacheckfile
USED="$(cat quotacheckfile | grep envhome -1 | tail -1 | gawk '{printf ("%.0f", $3)}')"
QUOTA="$(cat quotacheckfile | grep envhome -1 | tail -1 | gawk '{printf ("%.0f", $5)}')"
LEFT="$(expr $QUOTA - $USED)"
echo "disk space = $USED/$QUOTA; $LEFT left"
if [ $LEFT -lt 10 ]
then
   echo "remaing quota < 10 Gig: stopping resubmission of job"
fi