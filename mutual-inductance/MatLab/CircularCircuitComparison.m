function [results]...
    = CircularCircuitComparison(separations,... % An array of separation vectors (x,y,z) in meters
                                drivenRadius,... % Driven circuit radius in meters
                                measuredRadius) % Measured circuit radius in meters
                            
    faraday = FaradayEmfForCircularCircuits;
    eField = InductiveFieldEmfForCircularCircuits;
    
    dIdt = 1000.0;

    results.faradayEmfCoplanar = [];
    results.proposedAccelEmfCoplanar = [];
    results.conventionalAccelEmfCoplanar = [];
    
    separationForCoAxial = -(drivenRadius + measuredRadius);
    
    % Calculate coplanar orientations, so vertical separation is zero
    for i = 1:length(separations)
        sep = separations(i);
        emfFaraday = faradayEmf(faraday, dIdt, drivenRadius, measuredRadius, sep, 0);
        results.faradayEmfCoplanar(i) = -emfFaraday;
        
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;
        emfProposed = totalEmf(eField, dIdt, drivenRadius, measuredRadius, sep, 0);
        results.proposedAccelEmfCoplanar(i) = emfProposed;

        eField.proposedAccelTermCoeff = 0;
        eField.conventionalAccelTermCoeff = 1;
        emfConventional = totalEmf(eField, dIdt, drivenRadius, measuredRadius, sep, 0);
        results.conventionalAccelEmfCoplanar(i) = emfConventional;

    end

    results.faradayEmfCoaxial = [];
    results.proposedAccelEmfCoaxial = [];
    results.conventionalAccelEmfCoaxial = [];
    
    % Calculate coplanar orientations, so vertical separation is zero
    for i = 1:length(separations)
        sep = separations(i);
        emfFaraday = faradayEmf(faraday, dIdt, drivenRadius, measuredRadius, separationForCoAxial, sep);
        results.faradayEmfCoaxial(i) = -emfFaraday;
        
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;
        emfProposed = totalEmf(eField, dIdt, drivenRadius, measuredRadius, separationForCoAxial, sep);
        results.proposedAccelEmfCoaxial(i) = emfProposed;

        eField.proposedAccelTermCoeff = 0;
        eField.conventionalAccelTermCoeff = 1;
        emfConventional = totalEmf(eField, dIdt, drivenRadius, measuredRadius, separationForCoAxial, sep);
        results.conventionalAccelEmfCoaxial(i) = emfConventional;

    end
    
    m2mm = 1000;
    micro = 1.0e6;
    
    fprintf("                           Coplanar(uV)                                 Coaxial(uV)\n");
    fprintf("separation (mm)  Faraday    Prop accel   Conv accel           Faraday    Prop accel   Conv accel\n");
    for i = 1:length(separations)
    fprintf("% -10.4g       % -10.4g  % -10.4g  % -10.4g           % -10.4g  % -10.4g  % -10.4g\n", m2mm * separations(i),...
        micro * results.faradayEmfCoplanar(i),...
        micro * results.proposedAccelEmfCoplanar(i),...
        micro * results.conventionalAccelEmfCoplanar(i),...
        micro * results.faradayEmfCoaxial(i),...
        micro * results.proposedAccelEmfCoaxial(i),...
        micro * results.conventionalAccelEmfCoaxial(i));
    end
    
end
