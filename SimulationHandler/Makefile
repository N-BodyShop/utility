STRUCTURES_PATH = $(HOME)/projects/structures
CC = $(CXX)
CXXFLAGS += -g -O3 -Wall -I$(STRUCTURES_PATH) -I$(XERCESCROOT)/include
LDFLAGS += -L$(XERCESCROOT)/lib -L$(STRUCTURES_PATH) 
LDLIBS += -lxerces-c -lTipsy

TARGETS = simload

all:: $(TARGETS)

simload: SiXFormat.o TipsyFormat.o

TipsyFormat.o: Simulation.h

SiXFormat.o: Simulation.h

simload.o: Simulation.h SiXFormat.h TipsyFormat.h

clean::
	rm -f *~ *.o core.* $(TARGETS)
	
