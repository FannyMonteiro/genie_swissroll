#
# $Id: Makefile 3032 2006-11-27 18:44:58Z sp1003 $
#
##################################################################
#               CALLEE SIDE STUBS & LIBRARIES                    #
##################################################################

include ../../../genie-main/makefile.arc
include ../../../babel-main/makefile.arc

##################################################################
#                            SETTINGS                            #
##################################################################

LEGACY_IMPLEMENTATION_LIB  = ${GENIE_ROOT}/genie-embm/lib
LEGACY_MODULE              = embm

CASTPATH                   = ${GENIE_ROOT}/babel-main/lib

##################################################################
#                       DO NOT MODIFY                            #
##################################################################

INCL_GENIE_MAIN            = -I${GENIE_ROOT}/genie-main
CAST_FILES                 = castings.F90

include babel.make
include glue/babel.make

EXT_OBJS                   = ${BASICMODULESRC} \
                             ${TYPEMODULESRCS} \
                             ${ARRAYMODULESRCS} \
                             ${STUBMODULESRCS} \
                             ${IORSRCS:.c=.o} \
                             ${SKELSRCS:.c=.o} \
                             ${STUBSRCS:.c=.o}

make_dynamic_lib: gluefiles
	cd glue ; $(F90LINK) ${EXT_OBJS} ${CASTPATH}/${CAST_FILES} \
        ../${IMPLMODULESRCS} ../${IMPLSRCS} \
        -o ../libimplementation.so -L${LEGACY_IMPLEMENTATION_LIB} -l${LEGACY_MODULE} \
        ${GENIE_LIB_INCLUDE} ${INCL_GENIE_MAIN} ${NETCDF_INCLUDE} \
        ${LINK_LIBS_IMPL}

gluefiles:
	cd glue ; ${MAKE}

clean:
	rm -f *.o *.mod *.so
	cd glue ; ${MAKE} clean
