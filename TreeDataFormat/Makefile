STRUCTURES_PATH = $(HOME)/projects/structures
CC = $(CXX)
CXXFLAGS += -Wall -g -O3  -I$(STRUCTURES_PATH) #-I$(HOME)/projects/boost_1_29_0
LDFLAGS += -L$(STRUCTURES_PATH)
LDLIBS += -lTipsy

TARGETS = tipsy2tree see_header browse_field #check_tree region_query #convert_ascii_field 

all:: $(TARGETS)

tipsy2tree.o: SFC.h $(STRUCTURES_PATH)/libTipsy.a

tipsy2tree.o: Makefile

clean::
	rm -f *~ *.o core.* $(TARGETS)
	
