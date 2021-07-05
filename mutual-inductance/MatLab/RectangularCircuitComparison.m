function [results]...
    = RectangularCircuitComparison(separations,... % An array of separation vectors (x,y,z) in meters
                                   drvDimension,...
                                   msrDimension)
                            
    faraday = FaradayEmfForRectangularCircuits;
    eField = InductiveFieldEmfForRectangularCircuits;
    
    dIdt = 1000.0;

    results.faradayEmfCoplanar = [];
    results.proposedAccelEmfCoplanar = [];
    results.conventionalAccelEmfCoplanar = [];
    
    srcWires = RectangularCircuitsCommon.wiresFromDimensions(drvDimension, [0,0,0]);
    msrOffsetX = (drvDimension(1) + msrDimension(1))/2;

    % Calculate coplanar orientations, so vertical separation is zero
    for i = 1:length(separations)
        sep = separations(i);        
        msrWires = RectangularCircuitsCommon.wiresFromDimensions(msrDimension, [msrOffsetX + sep, 0, 0]);
        
        emfFaraday = faradaysEmfWireSetToWireSet(faraday, dIdt, srcWires, msrWires);
        results.faradayEmfCoplanar(i) = emfFaraday.emfTotal;
        
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;
        emfProposed = inlineEmfForWireSetsOnWireSets(eField, dIdt, srcWires, msrWires);
        results.proposedAccelEmfCoplanar(i) = emfProposed.emfTotal;

        eField.proposedAccelTermCoeff = 0;
        eField.conventionalAccelTermCoeff = 1;
        emfConventional = inlineEmfForWireSetsOnWireSets(eField, dIdt, srcWires, msrWires);
        results.conventionalAccelEmfCoplanar(i) = emfConventional.emfTotal;

    end

    results.faradayEmfCoaxial = [];
    results.proposedAccelEmfCoaxial = [];
    results.conventionalAccelEmfCoaxial = [];
    
    % Calculate coplanar orientations, so vertical separation is zero
    for i = 1:length(separations)
        sep = separations(i);
        
        msrWires =  RectangularCircuitsCommon.wiresFromDimensions(msrDimension, [0, 0, sep]);
        emfFaraday = faradaysEmfWireSetToWireSet(faraday, dIdt, srcWires, msrWires);
        results.faradayEmfCoaxial(i) = emfFaraday.emfTotal;
  
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;
        emfProposed = inlineEmfForWireSetsOnWireSets(eField, dIdt, srcWires, msrWires);
        results.proposedAccelEmfCoaxial(i) = emfProposed.emfTotal;

        eField.proposedAccelTermCoeff = 0;
        eField.conventionalAccelTermCoeff = 1;
        emfConventional = inlineEmfForWireSetsOnWireSets(eField, dIdt, srcWires, msrWires);
        results.conventionalAccelEmfCoaxial(i) = emfConventional.emfTotal;

    end
    
    m2mm = 1000;
    micro = 1.0e6;
    
    fprintf("                           Coplanar(uV)                                 Coaxial(uV)\n");
    fprintf("separation (mm)  Faraday    Prop accel   Conv accel           Faraday    Prop accel   Conv accel\n");
    for i = 1:length(separations)
    fprintf("%-10g       %-10g  %-10g %-10g         %-10g  %-10g  %-10g\n", m2mm * separations(i),...
        micro * results.faradayEmfCoplanar(i),...
        micro * results.proposedAccelEmfCoplanar(i),...
        micro * results.conventionalAccelEmfCoplanar(i),...
        micro * results.faradayEmfCoaxial(i),...
        micro * results.proposedAccelEmfCoaxial(i),...
        micro * results.conventionalAccelEmfCoaxial(i));
    end

end
