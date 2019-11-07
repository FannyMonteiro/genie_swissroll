#
# $Id: Makefile 3048 2006-11-27 19:14:40Z sp1003 $
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