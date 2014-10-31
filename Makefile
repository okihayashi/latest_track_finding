TARGETS=TrackFinding
SRCS = TrackFinding.cc LayerInf140328.cc Wirepos0.cc WireposEP_kai.cc WireposEP2_kai.cc WireposReverse.cc NeuralNet.cc Density_Cut.cc Distance.cc DrawDetector.cc
OBJS = $(SRCS:.cc=.o)

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

CXXFLAGS  = -Wall -O2 $(ROOTFLAGS)
CXXLIBS   = $(ROOTLIBS)

all: $(TARGETS)

$(TARGETS): $(OBJS)
	g++ -o $@ $^ $(CXXLIBS)

.cc.o:
	g++ -c $(CXXFLAGS) $<

.PHONY: clean
clean:	
	rm -rf $(OBJS) $(TARGETS)
