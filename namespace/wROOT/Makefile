CC=g++
INC_DIR= .
CFLAGS=-O -g -Wall -I$(INC_DIR) $(shell root-config --cflags)
DEPS= IceRayTracing.hh
SOURCES= IceRayTracing.cc uzairRayTrace.cc
OBJ= $(SOURCES:.cc = .o)
LDFLAGS=-L${GSLDIR}/lib $(shell root-config --glibs)
CXXFLAGS=-I${GSLDIR}/include
LDLIBS=-lgsl -lgslcblas $(shell root-config --libs)


uzairRayTrace: $(OBJ)
	$(CC) $(CFLAGS) -o uzairRayTrace $(OBJ) $(LDLIBS)
