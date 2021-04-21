# Inductive Field Supplemental Materials

For the paper *An inductive field component of the electric field that corresponds to Faraday's law*, the self inductance of a straight wire is determined numerically. The derivations are included in the supplimental text for the paper. The calculations were done in C++ and the source code is available here.

## Building the C++ code

The C++ programs can be built using the CMake build system. For users unfamiliar with building using CMake, instructions can be found [here](https://cmake.org/runningcmake/).

There is only one program complied as part of a full build, `wire-self-inductance`, or on Windows systems, `wire-self-inductance.exe`.

## Program arguments

The typical command line for calling this function is as follows, the values being those for the 1mm diameter wire case:
```
wire-self-inductance --wire-length-mm=2.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
```

where:
* wire-length-mm is the length of the wire segemnt in millimeters;
* wire-diameter-mm is the diameer of the wire in millimeters;
* num-mesh-rings is how many rings of triangular elements are created over the cross-section of the wire; and
* layer-thickness-mm is the thickness of each element in the direction parallel to the axis of the wire.

Note that layer-thickness-mm is the target thickness. If this does not divide evenly along the length of the wire a value will be used so that all elements are the same thickness.

## Program output
