# Makefile from Jianglai

ObjSuf        = o
SrcSuf        = C
DllSuf        = so
OutPutOpt     = -o  
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../include

# Work with Linux with egcs	
CXX           = g++ 
CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += -lSpectrum -lMinuit

OBJECTS = BetaSpectrum.o plotFierz.o
SOURCE = plotFierz

.PHONY: all
all: $(SOURCE)

$(SOURCE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(SOURCE) $(OBJECTS) $(LIBS)

BetaSpectrum.o: BetaSpectrum.hh
plotFierz.o: comparehist.hh

# -------------------------------------------------------------------------------
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
#%: %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
#	@echo "$@ done"

clean:
		@rm -f *.o *~  core $(SOURCE) *.pdf