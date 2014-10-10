TARGETS=rearrange
SRCS = rearrange.cc
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
