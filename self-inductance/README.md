# Straight wire self inductance

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

The output for the above case of 2mm wire with 1mm diameter is the three inductance values, namely Rosa, the proposed inductive fields inductance, and the conventional terms inductance. Some information about the mesh is also given.

```
Rosa inductance (nH): 0.625574
Proposed inductance (nH): 0.618806
Conventional inductance (nH): 0.319592
    Wire diameter (mm): 1
    Wire length (mm): 2
    Mesh rings: 6
    Number of polygons in x-y: 144
    Layer thickness requested (mm): 0.025
    Layer thickness used (mm): 0.025
```

## Command lines for cases in the paper

```
wire-self-inductance --wire-length-mm=2.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
wire-self-inductance --wire-length-mm=5.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
wire-self-inductance --wire-length-mm=10.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
wire-self-inductance --wire-length-mm=15.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
wire-self-inductance --wire-length-mm=20.0 --wire-diameter-mm=1.0 --num-mesh-rings=6 --layer-thickness-mm=0.025
```

## Notes about the numerical solution

The general approach to the derivation is covered in the paper. It is simply a the sum of the interactions between all of the elements.

The problem with the straight forward approach is that with 800 layers and 144 polygons per layer there are about 13 billion combinations. Given the constant diameter of the wire it is clear that the same pattern repeats itself. For example, the effect of layer 1 on layer 2 will be the same as 2 on 3. It is also clear that the sum of the a set of layers can be applied to another layer of with the same number of interactions.

The calculation done here is that a list of force for distance between layers is created. Then a vector of sums of layers is created. For example, one entry would have the layer 1 on itself. The next entry would be the sum of layer 1 on itself and layer 2, and so on.

To speed up the processing is broken up into threads, each thread calculating a set of layers.

The majority of the code involved is for mesh creation and the bookkeeping related to interlayer forces and how they relate to the final sum.

## Acknowledgements

This program uses the Eigen C++ library for vector math and the args.hxx library to simplify command line argument handling.

