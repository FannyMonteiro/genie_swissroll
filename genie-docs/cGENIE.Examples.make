#!/bin/bash
#
######################################
### SCIPT TO MAKE Examples PS FILE ###
######################################
#

latex cGENIE.Examples.tex
latex cGENIE.Examples.tex
dvips -o cGENIE.Examples.ps cGENIE.Examples.dvi
evince cGENIE.Examples.ps &
rm cGENIE.Examples.dvi
rm cGENIE.Examples.out
rm cGENIE.Examples.aux
rm cGENIE.Examples.toc
