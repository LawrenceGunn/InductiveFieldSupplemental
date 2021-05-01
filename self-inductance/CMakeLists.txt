cmake_minimum_required(VERSION 3.17)
project(wire-self-inductance)

set(CMAKE_CXX_STANDARD 17)

set(wire-self-inductance_hdr
        args.hxx
        CircularMesh.h
        RosaStraightWireSelfInductance.h
        MeshedStraightWireSelfInductance.h
        )

set(wire-self-inductance_src
        CircularMesh.cpp
        main.cpp
        MeshedStraightWireSelfInductance.cpp
        )

add_executable(wire-self-inductance
        ${wire-self-inductance_hdr}
        ${wire-self-inductance_src}
        )