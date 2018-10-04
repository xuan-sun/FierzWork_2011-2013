ObjSuf		= o
SrcSuf		= C
DllSuf		= so
OutPutOpt	= -o
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
ROOTCFLAGS	= $(shell root-config --cflags)

VPATH = ./:../include

CXX		= g++
CXXFLAGS	= -O2 -Wall -fPIC -std=c++11
LD		= g++
SOFLAGS		= -shared
LIBS		= $(ROOTLIBS) $(ROOTGLIBS)
CXXFLAGS       += $(ROOTCFLAGS) -I../include
LIBS	       += -lSpectrum -lMinuit

objects = Asymmetries.o SQLinterface.o SystematicCorrections.o MBUtils.o EvtRateHandler.o calibrationTools.o


.PHONY: all
all: MBAnalyzer.exe UK2SIMcomparison.exe

MBAnalyzer.exe: MBAnalyzer.o $(objects) BetaSpectrum.o
	$(CXX) $(CXXFLAGS) -o MBAnalyzer.exe MBAnalyzer.o $(objects) BetaSpectrum.o $(LIBS) 
	@echo "Done"

UK2MPMcomparison.exe: UK2MPMcomparison.o $(objects) 
	$(CXX) $(CXXFLAGS) -o UK2MPMcomparison.exe UK2MPMcomparison.o $(objects) $(LIBS) 
	@echo "Done"

UK2SIMcomparison.exe: UK2SIMcomparison.o $(objects) 
	$(CXX) $(CXXFLAGS) -o UK2SIMcomparison.exe UK2SIMcomparison.o $(objects) $(LIBS) 
	@echo "Done"

SuperSumComparison.exe: SuperSumComparison.o $(objects) 
	$(CXX) $(CXXFLAGS) -o SuperSumComparison.exe SuperSumComparison.o $(objects) $(LIBS) 
	@echo "Done"

calibrationTools.o: calibrationTools.hh

Asymmetries.o: Asymmetries.hh MBUtils.o

BetaSpectrum.o: BetaSpectrum.hh

EvtRateHandler.o: EvtRateHandler.hh SQLinterface.hh MBUtils.o

SQLinterface.o: SQLinterface.hh

SystematicCorrections.o: SystematicCorrections.hh

MBUtils.o: MBUtils.hh

clean:
		@rm -f *.o *~  core
