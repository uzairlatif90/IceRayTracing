CC=g++
INC_DIR= .
CFLAGS=-O -g -Wall -I$(INC_DIR)
DEPS= IceRayTracing.hh
SOURCES= IceRayTracing.cc uzairRayTrace.cc
OBJ= $(SOURCES:.cc = .o)
LDFLAGS=-L${GSLDIR}/lib
CXXFLAGS=-I${GSLDIR}/include
LDLIBS=-lgsl -lgslcblas

uzairRayTrace: $(OBJ)
	$(CC) $(CFLAGS) -o uzairRayTrace $(OBJ) $(LDLIBS)
