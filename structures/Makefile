CXXFLAGS += -g -O3 -Wall

NAME = structures
VERSION = 1.1
INSTALL_LOCATION = $(HOME)/my_stuff/bin/
LIBTARGETS = TipsyReader.o TipsyFile.o SS.o
EXECTARGETS = browse see_header random_sample

.PHONY: all docs dist clean

all: libTipsy.a($(LIBTARGETS)) SFC.o $(EXECTARGETS)

TipsyFile.o: TipsyFile.h TipsyParticles.h TipsyFile.cpp Vector3D.h OrientedBox.h

TipsyReader.o: TipsyReader.h TipsyReader.cpp TipsyParticles.h Vector3D.h xdr_template.h

SS.o: xdr_template.h Vector3D.h

SFC.o: SFC.cpp SFC.h
	$(CXX) -O -ffloat-store -c SFC.cpp

browse see_header random_sample: libTipsy.a

docs:
	doxygen structures.doxygen

dist:
	mkdir $(NAME)-$(VERSION)
	cp Makefile *.h *.cpp structures.doxygen $(NAME)-$(VERSION)/
	tar zcf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	rm -Rf $(NAME)-$(VERSION)
	
clean:
	rm -f core* *.o *~ $(TARGETS) libTipsy.a

install: all
	cp $(EXECTARGETS) $(INSTALL_LOCATION)
