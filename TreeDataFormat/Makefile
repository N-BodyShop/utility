STRUCTURES_PATH = $(HOME)/projects/structures
CC = $(CXX)
CXXFLAGS += -Wall -O3 -I$(STRUCTURES_PATH)
LDFLAGS += -L$(STRUCTURES_PATH)
LDLIBS += -lTipsy $(STRUCTURES_PATH)/SFC.o -lpopt

TARGETS = tipsy2tree see_header browse_field check_tree convert_to_ascii histogram_field

all:: $(TARGETS)

tipsy2tree.o: $(STRUCTURES_PATH)/libTipsy.a

tipsy2tree.o: Makefile

clean::
	rm -f *~ *.o core.* $(TARGETS)
	
