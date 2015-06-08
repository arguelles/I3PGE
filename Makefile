CXX = clang++
CXXFLAGS+=-O3 -std=c++11 -fPIC

PHOTOSPLINE=/home/carguelles/programs/photospline_hack/
SQUIDS=/home/carguelles/programs/SQuIDS/build/
nuSQUIDS=/home/carguelles/programs/nuSQuIDS/build/

I3PGE_PATH=${CURDIR}

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

CXXFLAGS+=-I$(PHOTOSPLINE)/src/include
CXXFLAGS+=-I$(nuSQUIDS)/include -I$(SQUIDS)/include
CXXFLAGS+=-I/data/user/cweaver/tools/RHEL_6_x86_64/include/
CXXFLAGS+=-I./inc

# Project files
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

#Dynamic Library
OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,lib$(NAME).so
endif
ifeq ($(OS_NAME),Darwin)
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -install_name $(I3PGE_PATH)/lib/$(DYN_PRODUCT) -compatibility_version $(VERSION) -current_version $(VERSION)
endif
VERSION=1.0.0
NAME=I3PGE
STAT_PRODUCT=lib$(NAME).a
DYN_PRODUCT=lib$(NAME)$(DYN_SUFFIX)

.PHONY: all clean

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

examples:
	cd examples; make;

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)
	mv $(DYN_PRODUCT) $(I3PGE_PATH)/lib/$(DYN_PRODUCT)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)
	mv $(STAT_PRODUCT) $(I3PGE_PATH)/lib/$(STAT_PRODUCT)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o lib/* bin/*

