CC=g++
INC_DIR= .
CFLAGS=-O -g -Wall -I$(INC_DIR) `root-config --cflags`
DEPS= IceRayTracing.hh
SOURCES= IceRayTracing.cc uzairRayTrace.cc
OBJ= $(SOURCES:.cc = .o)
LDFLAGS=-L${GSLDIR}/lib `root-config --glibs`
CXXFLAGS=-I${GSLDIR}/include
LDLIBS=-lgsl -lgslcblas `root-config --libs`


uzairRayTrace: $(OBJ)
	$(CC) $(CFLAGS) -o uzairRayTrace $(OBJ) $(LDLIBS)
