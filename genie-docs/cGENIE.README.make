#!/bin/bash
#
####################################
### SCIPT TO MAKE README PS FILE ###
####################################
#

latex cGENIE.README.tex
dvips -o cGENIE.README.ps cGENIE.README.dvi
evince cGENIE.README.ps &
rm cGENIE.README.dvi
rm cGENIE.README.out
rm cGENIE.README.aux
rm cGENIE.README.toc
