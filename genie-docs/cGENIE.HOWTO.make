#!/bin/bash
#
###################################
### SCIPT TO MAKE HOWTO PS FILE ###
###################################
#

latex cGENIE.HOWTO.tex
latex cGENIE.HOWTO.tex
dvips -o cGENIE.HOWTO.ps cGENIE.HOWTO.dvi
evince cGENIE.HOWTO.ps &
rm cGENIE.HOWTO.dvi
rm cGENIE.HOWTO.out
rm cGENIE.HOWTO.aux
rm cGENIE.HOWTO.toc
