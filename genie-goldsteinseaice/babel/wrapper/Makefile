#
# $Id: Makefile 3049 2006-11-27 19:16:40Z sp1003 $
#
##################################################################
#               CALLER SIDE STUBS & LIBRARIES                    #
##################################################################

include ../../../genie-main/makefile.arc
include ../../../babel-main/makefile.arc

CALLER_PATH        = ${GENIE_ROOT}/genie-seaice/babel/impl

##################################################################
#                          SETTINGS                              #
##################################################################

WRAPPER_CLASS_DIR  = seaice
WRAPPER_LIB        = libseaice.so

##################################################################
#                       DO NOT MODIFY                            #
##################################################################

BABEL_LIB_DIR      = ${BABEL_HOME}/lib/
INCLUDE_PATHS      = ${BABEL_LIB_DIR}/sidl-0.11.0.jar:${BABEL_LIB_DIR}/sidlstub_0.11.0.jar

include glue/babel.make

wrapper: make_dynamic_lib
	cd ${WRAPPER_CLASS_DIR} ; \
	${JAVA_HOME}/bin/javac -classpath ${CLASSPATH}:${INCLUDE_PATHS} *_interface.java ; \
	${JAVA_HOME}/bin/javac -classpath ${CLASSPATH}:${INCLUDE_PATHS} *.java

make_dynamic_lib:
	cd glue/ && ${LINKER} ${STUBSRCS} \
        -o ../${WRAPPER_LIB} ${CALLER_PATH}/libimplementation.so \
        `${BABEL_HOME}/bin/babel-config --libs-f90` \
        `${BABEL_HOME}/bin/babel-config --includes` \
        ${JNI_INCLUDES}

clean:
	rm -rf **/*.class *.so glue/*.o glue/*.so glue/lib*

