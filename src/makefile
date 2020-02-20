CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11

GRID_SRC = grid.cc grid_test.cc
GRID_INC = grid.h
GRID_OBJ = $(patsubst %.cc,%.o,$(GRID_SRC))

VTK_SRC = vtk_writer.cc
VTK_INC = vtk_writer.h
VTK_OBJ = $(patsubst %.cc,%.o,$(VTK_SRC))

INC = $(GRID_INC) $(VTK_INC)
OBJ = $(GRID_OBJ) $(VTK_OBJ)

all: main ccfv1

#------------------------------------------------------------------------------
grid.o: $(GRID_INC) grid.cc
	$(CXX) $(CXXFLAGS) -c grid.cc

grid_test.o: $(GRID_INC) grid_test.cc
	$(CXX) $(CXXFLAGS) -c grid_test.cc
#------------------------------------------------------------------------------
vtk_writer.o: $(GRID_INC) $(VTK_INC) vtk_writer.cc
	$(CXX) $(CXXFLAGS) -c vtk_writer.cc
#------------------------------------------------------------------------------
main.o: $(INC) main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

main: $(OBJ) main.o
	$(CXX) -o main main.o $(OBJ)
#------------------------------------------------------------------------------
ccfv1.o: $(INC) ccfv1.cc
	$(CXX) $(CXXFLAGS) -c ccfv1.cc

ccfv1: $(OBJ) ccfv1.o
	$(CXX) -o ccfv1 ccfv1.o $(OBJ)
#------------------------------------------------------------------------------

clean:
	rm -rf *.o main
