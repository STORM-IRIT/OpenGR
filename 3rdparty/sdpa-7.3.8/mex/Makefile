#
# SDPA-M: $Revision: 7.3 $
#
# For Matlab
# $  make
# For Octave
# $ make COMPILE_ENVIRONMENT=octave

# after "make install", you can find
# 'make.inc' in 'share/sdpa' sub-directory under the installed directory
# For example
# make MAKE_INCLUDE_DIR=/usr/share/sdpa
#  or 
# make MAKE_INCLUDE_DIR=/usr/local/share/sdpa
# In addition, for Octave, you should try
# make MAKE_INCLUDE_DIR=/usr/local/share/sdpa COMPILE_ENVIRONMENT=octave


# To compile with mingw, the following command is a sample.
# Your need to modify directory names and compiler names.
# make MEX="i586-mingw32msvc-g++ -shared" COMPILE_ENVIRONMENT=octave \
# PRINTF_INT_STYLE=-DPRINTF_INT_STYLE=\\\"%zd\\\" \
# OUTPUT_FORMAT="-o \$@.mexw32" \
# MATLAB_INCLUDE="-Imatlab-extern/include" \
# MATLAB_LIBS="-I/matlab-extern/microsoft/libmx.lib -I/matlab-extern/microsoft/libmex.lib" \


# for only SDPA-Project members
# set DEVELOPEMENT_STAGE=1 to compile without 'make install'
# For example
# make MAKE_INCLUDE_DIR=../etc DEVELOPEMENT_STAGE=1 
# In addition, for Octave, you should try
# make MAKE_INCLUDE_DIR=../etc DEVELOPEMENT_STAGE=1 COMPILE_ENVIRONMENT=octave

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

MAKE_INCLUDE_DIR=..
-include ${MAKE_INCLUDE_DIR}/make.inc
# after "make install", you can find
# 'make.inc' in 'share/sdpa' sub-directory under the installed directory
# For example
# make MAKE_INCLUDE_DIR=/usr/share/sdpa
#  or 
# make MAKE_INCLUDE_DIR=/usr/local/share/sdpa


# Set manually the following MEX full path
MEX = "auto"
# for Linux (depends on version)
#MEX = /usr/local/matlab/MatlabR2008b/bin/mex
# for Mac (depends on version)
OS = $(shell uname)
ifeq ($(OS),Darwin)
	MEX = /Applications/MATLAB_R2008b.app/bin/mex
endif

# ---- USUALLY just 'mex' is used
MEX = mex

# Automatic search MEX full path by invoking matlab
ifeq ($(MEX),"auto")
	MEX = ${shell matlab -nodisplay -r "s=matlabroot; fprintf('%s\n',s);exit" | tail -2 | head -1 }/bin/mex
endif

MATLAB_LAPACK = -lmwlapack -lmwblas
DEF_FPIC = -fPIC
ARCH = $(shell uname -m)

CXXFLAGSOPTIONS = CXXFLAGS="-Wall  ${DEF_FPIC}" -O -DNDEBUG CC=${CC} CXX=${CXX}
LARGEARRAYOPTION = -largeArrayDims

SDPA_INCLUDE = -I${SDPA_DIR}/include
ifeq ($(DEVELOPEMENT_STAGE),1)
	SDPA_INCLUDE = -I..
	SDPA_LIB = ../libsdpa.a
endif
ALL_INCLUDE = ${SDPA_INCLUDE} ${MUMPS_INCLUDE} ${PTHREAD_INCLUDE}
# ALL_LIBS    = ${SDPA_LIB} ${MUMPS_LIBS} ${MATLAB_LAPACK} ${FCLIBS}
ALL_LIBS    = ${SDPA_LIB} ${MUMPS_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${PTHREAD_LIBS} ${FCLIBS}

PRINTF_INT_STYLE = -DPRINTF_INT_STYLE=\\\"%zd\\\"
# Output format is usually automatically set by mex command.
OUTPUT_FORMAT =

ifeq ($(COMPILE_ENVIRONMENT),octave)
        MEX = CC=${CC} CXX=${CXX} mkoctfile --mex
        CXXFLAGSOPTIONS =
        LARGEARRAYOPTION =
	PRINTF_INT_STYLE = -DPRINTF_INT_STYLE=\\\"%d\\\"
endif

ALL_OPTION = ${CXXFLAGSOPTIONS} ${LARGEARRAYOPTION} ${PRINTF_INT_STYLE}

all: mexsdpa mexSedumiWrap mexAggSDPcones mexDisAggSDPsol mexWriteSedumiToSDPA mexReadSDPAToSedumi mexReadOutput
# all: mexSedumiWrap

mexBinaryCheck:
	@echo "---------------------------------------------"
	@echo "Mex Compilation Start                        "
	@echo "---------------------------------------------"
	@if [ -f $(MEX) ]; then \
	echo "---------------------------------------------" ;\
	echo " We use                                      " ;\
	echo "   " $(MEX) ;\
	echo "   for Mex compilaton                        " ;\
	echo "---------------------------------------------" ;\
	else \
	echo "---------------------------------------------" ;\
	echo "Automatic search is failed,                  " ;\
	echo " Set the full path for mex command mannually " ;\
	echo " in "`pwd`"/Makekefile" ;\
	echo " and try 'make' in `pwd`                     " ;\
	echo "---------------------------------------------" ;\
	false ;\
	fi

mexsdpa:
	${MEX} ${ALL_OPTION} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_INCLUDE} ${MATLAB_INCLUDE} \
	  mexsdpa.cpp mexFprintf.c \
	  ${MATLAB_LIBS} ${ALL_LIBS} 

mexSedumiWrap:
	${MEX} ${ALL_OPTION} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_INCLUDE}  ${MATLAB_INCLUDE} \
	  mexSedumiWrap.cpp mexFprintf.c \
	  ${MATLAB_LIBS} ${ALL_LIBS} 

mexAggSDPcones:
	${MEX} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_OPTION}  ${MATLAB_INCLUDE} \
	  mexAggSDPcones.cpp mexFprintf.c ${MATLAB_LIBS}

mexDisAggSDPsol:
	${MEX} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_OPTION}  ${MATLAB_INCLUDE} \
	  mexDisAggSDPsol.cpp mexFprintf.c ${MATLAB_LIBS}

mexWriteSedumiToSDPA:
	${MEX} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_OPTION}  ${MATLAB_INCLUDE} \
	  mexWriteSedumiToSDPA.cpp ${MATLAB_LIBS}
mexReadSDPAToSedumi:
	${MEX} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_OPTION}  ${MATLAB_INCLUDE} \
	  mexReadSDPAToSedumi.cpp ${MATLAB_LIBS}
mexReadOutput:
	${MEX} \
	  ${OUTPUT_FORMAT} \
	  ${ALL_OPTION}  ${MATLAB_INCLUDE} \
	  mexReadOutput.cpp ${MATLAB_LIBS}

# clean all generated files
clean:
	rm -f *~
	rm -f mexsdpa.mex* mexSedumiWrap.mex*
	rm -f mexAggSDPcones.mex* mexDisAggSDPsol.mex*
	rm -f mexFprintf.mex*
	rm -f mexWriteSedumiToSDPA.mex*
	rm -f mexReadSDPAToSedumi.mex* mexReadOutput.mex*
distclean: clean

#
# End of File
#
