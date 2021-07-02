classdef (ConstructOnLoad) InductiveFieldEmfForRectangularCircuits
    properties
        e0 = 8.85418782e-12
        c = 299792458
        proposedAccelTermCoeff = 1.0
        conventionalAccelTermCoeff = 1.0
        m0
    end
    methods
        function obj = InductiveFieldEmfForRectangularCircuits
            obj.m0 = 1.0/(obj.e0 * obj.c^2);
        end
        
        % This function calculates the numerator of the acceleration terms.
        % The choice of which to use, or to use both, is from the term
        % coefficients, which are typically 1 and 0 or 0 and 1,
        % respectively.
        function prod = accelProduct(obj, accel, r)
            rUnit = r/norm(r);
            prod = obj.conventionalAccelTermCoeff * cross(rUnit, cross(rUnit, accel))...
                - obj.proposedAccelTermCoeff * accel;
        end

        function emf = emfAtPoint(obj, dIdtVec, rDrv, rMsr)
            rMsrDrv = rMsr - rDrv;
            emf = accelProduct(obj, dIdtVec, rMsrDrv) / (4.0 * pi * obj.e0 * norm(rMsrDrv) * obj.c^2);
        end
        
        function obj = set.proposedAccelTermCoeff(obj, coeff)
            obj.proposedAccelTermCoeff = coeff;
        end
        
        function obj = set.conventionalAccelTermCoeff(obj, coeff)
            obj.conventionalAccelTermCoeff = coeff;
        end

    end
    methods (Static)
        function vec = rVector(rDrvOuter, ... % Driven circuit radius (m)
                               rMsrOuter,... % Measured circuit radius (m)
                               separation,... % Distance between nearest points of the wire (m)
                               vertSeparation,... % Vertical separation of circuits (m). 0 for coplanar circuits
                               drvCircuitAngle,... % The angle around the driven circuit in radians
                               msrCircuitAngle) % The angle around the measured circuit in radians

                           vec = [
                               rDrvOuter * (1.0 - cos(drvCircuitAngle)) + separation + rMsrOuter * (1.0 + cos(msrCircuitAngle)),...
                               rMsrOuter * sin(msrCircuitAngle) - rDrvOuter * sin(drvCircuitAngle),...
                               vertSeparation];
        end

    end
end
                       