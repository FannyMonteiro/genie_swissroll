##################################################################
#                     ALL STUBS & LIBRARIES                      #
##################################################################

include ../../babel-main/makefile.arc

SIDL_FILE = embm.sidl

##################################################################
#                       MORE SETTINGS                            #
##################################################################

all:
	${BABEL_HOME}/bin/babel -sF90 -oimpl -u ${SIDL_FILE}
	${BABEL_HOME}/bin/babel -cJava -owrapper -u ${SIDL_FILE}
	cd impl ; ${MAKE}
	cd wrapper ; ${MAKE}

clean:
	cd impl ; ${MAKE} clean
	cd wrapper ; ${MAKE} clean