function result = RectangularCircuitsPaperExample()
    separations = [0.005 0.010 0.025 0.075]; % separation, horizontal and then vertical, in meters
    drvWireDim = [0.1, 0.1];                 % width and height of driven circuit in meters
    msrWireDim = [0.05, 0.04];               % width and height of measured circuit in meters

    fprintf("\n");
    RectangularCircuitComparison(separations, drvWireDim, msrWireDim);
    
    fprintf("\n\n");

    result = RectangularCircuitWireEmfs(separations(1), drvWireDim, msrWireDim);

    fprintf("\n");
end