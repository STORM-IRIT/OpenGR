MUMPS_VER = 4.10.0.dfsg
MUMPS_TAR_FILE = mumps_${MUMPS_VER}.orig.tar.gz

SDPA_MAKEFILE = ../Makefile
# This file is automatically made by './configure' script
###################################################################
# The following variables are set from
# the arguments given to './configure' script for building SDPA
###################################################################

CC          = ${shell grep "^CC ="        ${SDPA_MAKEFILE} | cut -f2 -d=}
CFLAGS      = ${shell grep "^CFLAGS = "   ${SDPA_MAKEFILE} | cut -f2 -d=}
CXX         = ${shell grep "^CXX ="       ${SDPA_MAKEFILE} | cut -f2 -d=}
CXXFLAGS    = ${shell grep "^CXXFLAGS = " ${SDPA_MAKEFILE} | cut -f2 -d=}
FC          = ${shell grep "^FC ="        ${SDPA_MAKEFILE} | cut -f2 -d=}
FCFLAGS     = ${shell grep "^FCFLAGS = "  ${SDPA_MAKEFILE} | cut -f2 -d=}
PTHREAD_INCLUDE = ${shell grep "^PTHREAD_INCLUDE = "  ${SDPA_MAKEFILE} | cut -f2 -d=}
AR          = ${shell grep "^AR = "       ${SDPA_MAKEFILE} | cut -f2 -d=}
RANLIB      = ${shell grep "^RANLIB = "  ${SDPA_MAKEFILE} | cut -f2 -d=}

###################################################################

all: build/lib/libdmumps.a

build/lib/libdmumps.a:${MUMPS_TAR_FILE}
	rm -rf build;
	tar xzf ${MUMPS_TAR_FILE}
	mv mumps-${MUMPS_VER}/ build;
	cd build; cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc;

	cd build; echo "CC = "   ${CC}       >> Makefile.inc;
	cd build; echo "FC = "   ${FC}       >> Makefile.inc;
	cd build; echo "FL = "   ${FC}       >> Makefile.inc;
	cd build; echo "OPTC = " ${CFLAGS}  ${PTHREAD_INCLUDE} >> Makefile.inc;
	cd build; echo "OPTF = " ${FCFLAGS} ${PTHREAD_INCLUDE}  >> Makefile.inc;
	cd build; echo "OPTL = " ${FCFLAGS} ${PTHREAD_INCLUDE}  >> Makefile.inc;
	cd build; echo "LPORDDIR = "`pwd`"/PORD/lib/"     >> Makefile.inc;
	cd build; echo "IPORDDIR = "`pwd`"/PORD/include/" >> Makefile.inc;
	cd build; echo "AR = " ${AR} rv " " >> Makefile.inc;
	cd build; echo "RANLIB = " ${RANLIB} >> Makefile.inc;
	cd build ; ${MAKE} d

clean:
	@echo "##########################################################"
	@echo "If you want to rebuild mumps dir, use 'make distclean'"
	@echo "##########################################################"

${MUMPS_TAR_FILE}:
	@echo "By wget command, downloading a file from"
	@echo "http://ftp.de.debian.org/debian/pool/main/m/mumps/${MUMPS_TAR_FILE}"
	@echo "If fail, you manually put the file into " `pwd` "directory."
	wget http://ftp.de.debian.org/debian/pool/main/m/mumps/${MUMPS_TAR_FILE}

distclean: 
	rm -rf build ${MUMPS_TAR_FILE}

install: all
uninstall:

