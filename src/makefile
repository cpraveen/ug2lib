CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11
GRID_SRC = grid.cc grid_test.cc
GRID_INC = grid.h
GRID_OBJ = $(patsubst %.cc,%.o,$(GRID_SRC))

all: main

grid.o: $(GRID_INC) grid.cc
	$(CXX) $(CXXFLAGS) -c grid.cc

grid_test.o: $(GRID_INC) grid_test.cc
	$(CXX) $(CXXFLAGS) -c grid_test.cc

main.o: $(GRID_INC) main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

main: $(GRID_OBJ) main.o
	$(CXX) -o main main.o $(GRID_OBJ)

clean:
	rm -rf *.o main
