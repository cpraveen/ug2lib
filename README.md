# A library for handling 2-D unstructured grids

At present the library supports reading only gmsh grids consisting of linear elements.

* `src/main.cc `: test grid preprocessor
* `src/ccfv1.cc`: first order FV code
* `src/ccfv2.cc`: second order FV code

## Formatting using astyle

Run astyle on source code before committing it
```
astyle --options=./astyle.rc grid.cc
```

## Authors

 * Shashwat Tiwari
 * Praveen Chandrashekar
