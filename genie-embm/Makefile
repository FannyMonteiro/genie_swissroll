#
# $Id: Makefile 6160 2011-03-21 10:52:03Z andy-ridgwell $
#

export

LIBNAME = libembm
FOR_DIR = src/fortran
LIB_DIR = lib
DATA_DIR = data

.PHONY : all $(LIBNAME).$(LIB_EXT) clean tags data

all: $(LIBNAME).$(LIB_EXT)

$(LIBNAME).$(LIB_EXT) :
	\cd $(FOR_DIR) ; $(MAKE) all
	\cd $(LIB_DIR) ; $(MAKE) all

clean :
	\cd $(FOR_DIR) ; $(MAKE) clean
	\cd $(LIB_DIR) ; $(MAKE) clean

tags:
	\cd $(FOR_DIR) ; $(MAKE) tags

data:
	\cd $(DATA_DIR)/download/; $(MAKE)
	\cd $(DATA_DIR)/input/; $(MAKE)
