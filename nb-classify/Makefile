# Makefile for TaxonScore

BINDIR = ../bin
OBJDIR = ../obj

CXX = g++
CXXFLAGS = -Wall -O3 -march=core2 -mfpmath=sse -msse2

COMPILE = $(CXX) $(CXXFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: nb-classify

nb-classify: $(OBJFILES)
	$(CXX) -o nb-classify $(OBJFILES)

%.o: %.cpp 
	$(COMPILE) -o $@ $<

clean:
	rm -f nb-classify *.o