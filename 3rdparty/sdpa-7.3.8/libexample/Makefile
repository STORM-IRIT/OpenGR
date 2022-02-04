
#############################################################
#  List of examples for callable-library
#
#  example1.cpp: Solve example1.dat-s using inputElements
#  example2.cpp: Solve example2.dat   using inputElements
#  example3.cpp: Solve example1.dat-s with initial point
#  example4.cpp: Solve twice with smaller epsilonStar
#  example5.cpp: Solve problems reading from file
#  example6.cpp: Solve problems reading from file with initial point
#
#############################################################

.PHONY: all lib clean cleanall distclean
.SUFFIXES: .exe

# after "make install", you can find
# 'make.inc' in 'share/sdpa' sub-directory under the installed directory
# For example
# make MAKE_INCLUDE_DIR=/usr/share/sdpa
#  or 
# make MAKE_INCLUDE_DIR=/usr/local/share/sdpa

MAKE_INCLUDE_DIR=..
-include ${MAKE_INCLUDE_DIR}/make.inc
# after "make install", you can find
# 'make.inc' in 'share/sdpa' sub-directory under the installed directory
# For example
# make MAKE_INCLUDE_DIR=/usr/share/sdpa
#  or 
# make MAKE_INCLUDE_DIR=/usr/local/share/sdpa


SRC = $(wildcard example?.cpp)
EXE = $(subst .cpp,.exe,$(SRC))

all: ${EXE}

%.exe: %.o
	${CXX} ${CXXFLAGS} -o $@ $< ${SDPA_LIBS} 
.cpp.o:
	${CXX} -c ${CXXFLAGS} \
		-I${SDPA_DIR}/include ${MUMPS_INCLUDE} ${PTHREAD_INCLUDE} \
	        -o $@ $<
clean:
	rm -f *.o *~
cleanall: clean
	rm -f *.exe
distclean: cleanall

install: all
uninstall:
