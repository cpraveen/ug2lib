GEO = $(wildcard *.geo)
MSH = $(GEO:.geo=.msh)

all: $(MSH)

%.msh: %.geo
	gmsh -2 -format msh22 $<

clean:
	rm -f $(MSH)
