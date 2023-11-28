if HAVE_PACKAGE_MUMPS
  SUBDIRS = etc 
else
  SUBDIRS = etc mumps
endif

includedir     = ${prefix}/include
datadir        = ${prefix}/share/@PACKAGE_NAME@
etcdir         = ${prefix}/share/@PACKAGE_NAME@/etc
docdir         = ${prefix}/share/@PACKAGE_NAME@/doc
exampledir     = ${prefix}/share/@PACKAGE_NAME@/example
libexampledir  = ${prefix}/share/@PACKAGE_NAME@/libexample
mexdir         = ${prefix}/share/@PACKAGE_NAME@/mex

man_MANS = sdpa.1
data_DATA = param.sdpa etc/make.inc

include_HEADERS = sdpa_algebra.h sdpa_block.h \
sdpa_call.h sdpa_chordal.h sdpa_dataset.h sdpa_dpotrf.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h

lib_LIBRARIES = libsdpa.a
libsdpa_a_SOURCES = sdpa_block.cpp sdpa_call.cpp sdpa_chordal.cpp \
sdpa_dataset.cpp sdpa_dpotrf.cpp sdpa_io.cpp sdpa_jordan.cpp \
sdpa_linear.cpp sdpa_newton.cpp sdpa_parts.cpp sdpa_solve.cpp \
sdpa_struct.cpp sdpa_tool.cpp \
sdpa_algebra.h sdpa_block.h \
sdpa_call.h sdpa_chordal.h sdpa_dataset.h sdpa_dpotrf.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h  
libsdpa_a_CXXFLAGS =  $(pthread_cflags) $(MUMPS_INCLUDE) $(PTHREAD_INCLUDE)

bin_PROGRAMS = sdpa
sdpa_SOURCES = sdpa_exe.cpp
sdpa_LDADD =  -L. -lsdpa $(MUMPS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(PTHREAD_LIBS) $(FCLIBS)
sdpa_CXXFLAGS = $(pthread_cflags) $(MUMPS_INCLUDE) $(PTHREAD_INCLUDE)

# each file of mumps is copied by 'nobase' option
#   due to sub-directory
if HAVE_PACKAGE_MUMPS
  # nothing to install 
else # install mumps
mumpsdir     = ${prefix}/share/@PACKAGE_NAME@

nobase_mumps_HEADERS = mumps/build/include/dmumps_c.h \
mumps/build/include/dmumps_root.h \
mumps/build/include/dmumps_struc.h \
mumps/build/include/mumps_c_types.h \
mumps/build/include/mumps_compat.h

nobase_mumps_DATA = mumps/build/lib/libdmumps.a \
mumps/build/lib/libmumps_common.a \
mumps/build/lib/libpord.a \
mumps/build/libseq/libmpiseq.a

endif # end of 'if HAVE_PACKAGE_MUMPS'

example_DATA = example1.dat-s example1.dat example1.ini \
example1.ini-s example2.dat

libexample_DATA  = libexample/Makefile \
libexample/example1.cpp libexample/example2.cpp \
libexample/example3.cpp libexample/example4.cpp \
libexample/example5.cpp libexample/example6.cpp

mex_DATA = mex/Makefile mex/SDPAToSedumi.m mex/SedumiToSDPA.m \
mex/gensdpafile.m mex/initial_point.m mex/mexAggSDPcones.cpp \
mex/mexDisAggSDPsol.cpp mex/mexFprintf.c mex/mexReadOutput.cpp \
mex/mexReadSDPAToSedumi.cpp mex/mexSedumiWrap.cpp \
mex/mexWriteSedumiToSDPA.cpp mex/mexsdpa.cpp mex/param.m \
mex/read_data.m mex/read_output.m mex/read_outputSedumi.m \
mex/sdpam.m mex/sedumiwrap.m \
mex/CommandList.txt

MOSTLYCLEANFILES = mex/*.mex* mex/*.o
