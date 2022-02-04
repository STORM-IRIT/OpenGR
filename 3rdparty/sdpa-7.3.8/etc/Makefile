# This file is a component of SDPA
# Copyright (C) 2004-2012 SDPA Project

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

SDPA_MAKEFILE = ../Makefile
# This file is automatically made by './configure' script
###################################################################
# The following variables are set from 
# the arguments given to './configure' script for building SDPA
###################################################################
SDPA_PACKAGE= ${shell grep "^PACKAGE_STRING =" ${SDPA_MAKEFILE} | awk '{printf("%s.%s",$$3,$$4)}'}

SDPA_DIR    = ${shell grep "^prefix ="         ${SDPA_MAKEFILE} | cut -f2 -d= | sed 's+ ++'}
SDPA_LIB    = ${SDPA_DIR}/lib/libsdpa.a 
MUMPS_DIR   = ${shell grep "^MUMPS_DIR ="      ${SDPA_MAKEFILE} | cut -f2 -d= | sed 's+ ++'}
MUMPS_INCLUDE   = ${shell grep "^MUMPS_INCLUDE ="  ${SDPA_MAKEFILE} | cut -f2 -d= | sed 's+ ++'}
MUMPS_LIBS  = ${shell grep "^MUMPS_LIBS ="     ${SDPA_MAKEFILE} | cut -f2 -d= | sed 's+ ++'}
CC          = ${shell grep "^CC ="           ${SDPA_MAKEFILE} | cut -f2 -d=}
CFLAGS      = ${shell grep "^CFLAGS = "      ${SDPA_MAKEFILE} | cut -f2 -d=}
CXX         = ${shell grep "^CXX ="          ${SDPA_MAKEFILE} | cut -f2 -d=}
CXXFLAGS    = ${shell grep "^CXXFLAGS = "    ${SDPA_MAKEFILE} | cut -f2 -d=}
FC          = ${shell grep "^FC ="           ${SDPA_MAKEFILE} | cut -f2 -d=}
FCFLAGS     = ${shell grep "^FCFLAGS = "     ${SDPA_MAKEFILE} | cut -f2 -d=}
LAPACK_LIBS = ${shell grep "^LAPACK_LIBS = " ${SDPA_MAKEFILE} | cut -f2 -d=}
BLAS_LIBS   = ${shell grep "^BLAS_LIBS = "   ${SDPA_MAKEFILE} | cut -f2 -d=}
PTHREAD_INCLUDE= ${shell grep "^PTHREAD_INCLUDE = " ${SDPA_MAKEFILE} | cut -f2 -d=}
PTHREAD_LIBS= ${shell grep "^PTHREAD_LIBS = " ${SDPA_MAKEFILE} | cut -f2 -d=}
FCLIBS      = ${shell grep "^FCLIBS = "      ${SDPA_MAKEFILE} | cut -f2 -d=}
VERSION     = ${shell grep "^VERSION = "     ${SDPA_MAKEFILE} | cut -f2 -d=}
###################################################################

TARGET = make.inc
all: ${TARGET}
${TARGET}:
	rm -f make.inc
	echo "VERSION     = " ${VERSION}   >> ${TARGET}
	echo "SDPA_DIR    = " ${SDPA_DIR}  >> ${TARGET}
	echo "SDPA_LIB    = " ${SDPA_LIB}  >> ${TARGET} 
	echo "MUMPS_DIR   = " ${MUMPS_DIR} >> ${TARGET}
	echo "CC          = " ${CC} >> ${TARGET}
	echo "CFLAGS      = " ${CFLAGS} >> ${TARGET}
	echo "CXX         = " ${CXX} >> ${TARGET}
	echo "CXXFLAGS    = " ${CXXFLAGS} >> ${TARGET}
	echo "FC          = " ${FC} >> ${TARGET}
	echo "FCFLAGS     = " ${FCFLAGS} >> ${TARGET}
	echo "LAPACK_LIBS = " ${LAPACK_LIBS} >> ${TARGET}
	echo "BLAS_LIBS   = " ${BLAS_LIBS} >> ${TARGET}
	echo "PTHREAD_INCLUDE= " ${PTHREAD_INCLUDE} >> ${TARGET}
	echo "PTHREAD_LIBS= " ${PTHREAD_LIBS} >> ${TARGET}
	echo "FCLIBS      = " ${FCLIBS} >> ${TARGET}
	echo "MUMPS_INCLUDE = " ${MUMPS_INCLUDE} >> ${TARGET}
	echo "MUMPS_LIBS = "  ${MUMPS_LIBS} >> ${TARGET}
	echo "SDPA_LIBS = " ${SDPA_LIB} ${MUMPS_LIBS}\
	       ${LAPACK_LIBS} ${BLAS_LIBS} \
	       ${FCLIBS} ${PTHREAD_LIBS}  >> ${TARGET}

force: clean all
clean:
	rm -f make.inc *~
distclean: clean
install: all
uninstall:

