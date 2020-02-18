CXX=g++
CXXFLAGS=-O3 -Wall -std=c++11

all: main

grid.o: grid.h grid.cc
	$(CXX) $(CXXFLAGS) -c grid.cc

grid_test.o: grid.h grid_test.cc
	$(CXX) $(CXXFLAGS) -c grid_test.cc

main.o: grid.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

main: grid.o grid_test.o main.o
	$(CXX) -o main main.o grid.o grid_test.o

clean:
	rm -rf *.o main
