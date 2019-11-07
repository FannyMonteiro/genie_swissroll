#!/bin/bash
#
#############################################
### SCIPT TO MAKE QuickStartGuide PS FILE ###
#############################################
#

latex cGENIE.QuickStartGuide.tex
dvips -o cGENIE.QuickStartGuide.ps cGENIE.QuickStartGuide.dvi
evince cGENIE.QuickStartGuide.ps &
rm cGENIE.QuickStartGuide.dvi
rm cGENIE.QuickStartGuide.out
rm cGENIE.QuickStartGuide.aux
rm cGENIE.QuickStartGuide.toc
