CXX=g++  -std=c++11
INCDIR=.
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

all: CustomFit.o test

CustomFit.o: CustomFit.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c CustomFit.C

test: test.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -o $@ CustomFit.o test.C

clean:
	rm *.o test
