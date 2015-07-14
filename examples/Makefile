# Compiler
CC=clang
CXX=clang++

PHOTOSPLINE=/home/carguelles/programs/photospline_hack/
SQUIDS=/home/carguelles/programs/SQuIDS/build/
nuSQUIDS=/home/carguelles/programs/nuSQuIDS/build/

I3PGE_PATH=${CURDIR}/../

# boost flags
LDFLAGS+=-lboost_system-mt -lboost_iostreams-mt -lboost_filesystem-mt -lboost_regex-mt
# hdf5 flags
LDFLAGS+=-L/data/user/cweaver/tools/RHEL_6_x86_64/lib -lhdf5 -lhdf5_hl
LDFLAGS+=-lhdf5 -lhdf5_hl
# photospline flags
LDFLAGS+=-L$(PHOTOSPLINE)/build
LDFLAGS+=-lphotospline
# nusquids flags
LDFLAGS+=-L$(SQUIDS)/lib -L$(nuSQUIDS)/lib -lnuSQuIDS -lSQuIDS

CXXFLAGS= -std=c++11
CXXFLAGS+=-I$(PHOTOSPLINE)/src/include
CXXFLAGS+=-I$(nuSQUIDS)/include -I$(SQUIDS)/include
CXXFLAGS+=-I/data/user/cweaver/tools/RHEL_6_x86_64/include/
CXXFLAGS+=-I../inc

EXAMPLES_SRC=$(wildcard *.cpp)
EXAMPLES=$(EXAMPLES_SRC:.cpp=.exe)

examples : $(EXAMPLES)

%.exe : %.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS) $< $(LDFLAGS) -o $@
	mv $@ bin/

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.exe lib/* bin/*

