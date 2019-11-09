ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./	
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .C
.PREFIXES: ./obj/

#----------------------------------------------------#

all: convertUSBWC2root8ChannelsBin  convertUSBWCMeas2root8ChannelsBin  runL1Bin

runL1Bin: runL1Bin.cpp obj/wfSim.o obj/waveformSimple.o obj/anaL1Bin.o
	g++ -o runL1Bin runL1Bin.cpp $(OUTLIB)*.o `root-config --cflags --glibs`

obj/anaL1Bin.o: src/anaL1Bin.cpp src/anaL1Bin.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)anaL1Bin.o $<

obj/anaMeasures.o: src/anaMeasures.cpp src/anaMeasures.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)anaMeasures.o $<

obj/wfSim.o: wfSimulation/src/wfSim.cpp wfSimulation/src/wfSim.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)wfSim.o $<

obj/waveformSimple.o: wfSimulation/src/waveformSimple.cpp wfSimulation/src/waveformSimple.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)waveformSimple.o $<

obj/waveform.o: wfSimulation/src/waveform.cpp wfSimulation/src/waveform.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)waveform.o $<

convertUSBWC2root8ChannelsBin: convertUSBWC2root8ChannelsBin.cc
	g++ -o convertUSBWC2root8ChannelsBin convertUSBWC2root8ChannelsBin.cc `root-config --cflags --glibs`

convertUSBWCMeas2root8ChannelsBin: convertUSBWCMeas2root8ChannelsBin.cc
	g++ -o convertUSBWCMeas2root8ChannelsBin convertUSBWCMeas2root8ChannelsBin.cc `root-config --cflags --glibs`

clean:
	rm -f runL1
	rm -f runL1Bin
	rm -f convertUSBWC2root8ChannelsBin
	rm -f convertUSBWCMeas2root8ChannelsBin
	rm -f *~
	rm -f src/*~
	rm -f $(OUTLIB)*.o
