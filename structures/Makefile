CXXFLAGS = -O3 -Wall

NAME = structures
VERSION = 1.1
TARGETS = TipsyFile.o TipsyFile.debug.o
all: $(TARGETS) tools

TipsyFile.o: TipsyFile.h TipsyParticles.h TipsyFile.cpp Vector3D.h OrientedBox.h

TipsyFile.debug.o: TipsyFile.h TipsyFile.cpp Vector3D.h OrientedBox.h
	$(CXX) -Wall -g -o TipsyFile.debug.o -c TipsyFile.cpp

tools::
	cd tools; make

docs:
	doxygen structures.doxygen

dist:
	mkdir $(NAME)-$(VERSION)
	cp Makefile *.h TipsyFile.cpp structures.doxygen $(NAME)-$(VERSION)/
	tar zcf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	rm -Rf $(NAME)-$(VERSION)
	
clean:
	rm -f core* *.o *~ $(TARGETS)
