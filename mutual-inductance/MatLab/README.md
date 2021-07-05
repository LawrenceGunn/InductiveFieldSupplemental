# Mutual inductance of wire circuits: MATLAB scripts

For the paper *An inductive field component of the electric field that corresponds to Faraday's law*, the mutual inductance between two wire circuits is determined numerically using MATLAB. These correspond to the Mathematica notebooks which were used for calculating the values in the paper.

To execute the examples run the two scripts:

    RectangularCircuitsPaperExample;
    CircularCircuitsPaperExample;
    
The output will correspond to the first three tables of the paper. The values should differ slightly from those of Mathematica as they are numerical results. The default accuracy specifications were used for both Mathematica and MATLAB, the latter having a slightly higher default tolerance.

The scripts were developed in MATLAB version 2021a.
