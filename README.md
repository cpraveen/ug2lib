# A library for handling 2-D unstructured grids

At present the library supports reading only gmsh grids consisting of linear elements of triangles and quadrilaterals. The mesh can consist of both types of cells. It uses the older version of Gmsh grid format, so set

```
Mesh.MshFileVersion = 2.2;
```

in your geo file. Or specify the format on the command line

```
gmsh -2 -format msh22 foo.geo -o foo.msh
```

* `src/main.cc `: test grid preprocessor
* `src/ccfv1.cc`: first order FV code
* `src/ccfv2.cc`: second order FV code

## Test the grid preprocessor

```
cd src
make
cd ../tests
gmsh -2 cylinder.geo
../src/main
```

See the output files `cylinder.vtk` and `test_msh.msh`.

```
visit -o cylinder.vtk
gmsh test_msh.msh
```

## Test first order FV code

```
cd tests
gmsh -2 square_tri.geo -o ccfv1.msh
../src/ccfv1
visit -o sol*.vtk
```

## Test second order FV code

This reads some parameters from an input file, see `tests/input.txt` file.

## Formatting using astyle

Run astyle on source code before committing it
```
astyle --options=./astyle.rc grid.cc
```

## Authors

 * Shashwat Tiwari
 * Praveen Chandrashekar
