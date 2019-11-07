#
# $Id: Makefile 3034 2006-11-27 18:52:20Z sp1003 $
#
##################################################################
#                        GLUE STUBS                              #
#                       DO NOT MODIFY                            #
##################################################################

include ../../../../genie-main/makefile.arc

.SUFFIXES: .c .o

.c.o:
	${CC} `${BABEL_HOME}/bin/babel-config --includes` -c $<

include babel.make

EXT_OBJS            =  ${BASICMODULESRC} \
                       ${TYPEMODULESRCS} \
                       ${ARRAYMODULESRCS} \
                       ${STUBMODULESRCS} \
                       ${IORSRCS:.c=.o} \
                       ${SKELSRCS:.c=.o} \
                       ${STUBSRCS:.c=.o}

all: ${EXT_OBJS}

clean:
	rm -f *.o *.mod *.so