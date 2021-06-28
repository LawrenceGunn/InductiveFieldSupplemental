classdef (ConstructOnLoad) InductiveFieldEmfForCircularCircuits
    properties
        e0 = 8.85418782e-12
        c = 299792458
        proposedAccelTermCoeff = 1.0
        conventionalAccelTermCoeff = 1.0
        m0
    end
    methods
        function obj = InductiveFieldEmfForCircularCircuits
            obj.m0 = 1.0/(obj.e0 * obj.c^2);
        end
           
        function emf = emfAtPoint(obj,...
                                  dIdt,...
                                  rMsr,...
                                  rDrvOuter,...
                                  rMsrOuter,...
                                  separation,...
                                  vertSeparation,...
                                  drvCircuitAngle,...
                                  msrCircuitAngle)
                            
            rVec = FaradayEmfForCircularCircuits.rVector(rDrvOuter, rMsrOuter, separation, vertSeparation, rMsr, drvCircuitAngle, msrCircuitAngle);
            dIdtVec = dIdt * [-sin(drvCircuitAngle), cos(drvCircuitAngle), 0];
            emf = obj.m0 / (4 * pi) * rMsr * rDrvOuter * dot(cross(dIdtVec, rVec), [0, 0, 1]) / norm(rVec)^3;
        end
        
        % Calculate the EMF separated by a vector r. Use of the proposed 
        % or conventional term is based on the proposedAccelTermCoeff or
        % conventionalAccelTermCoeff set by the two setter functions below.
        function emf = emfBetweenPoints(obj, dIdtVec, r)
            rLength = norm(r);
            rUnit = r./rLength;
            convComp = obj.conventionalAccelTermCoeff * cross(rUnit, cross(rUnit, dIdtVec));
            propComp = -obj.proposedAccelTermCoeff * dIdtVec;
            emf = (convComp + propComp)/(4 * pi * obj.e0 * rLength * obj.c^2);
        end
           
        % This function calulates the total EMF in measured wire using the LW theory.
        function emf = totalEmf(obj,...
                                dIdt,... % Change of current A/s
                                rDrvOuter,... %Driven circuit radius (m)
                                rMsrOuter,... % Measured circuit radius (m)
                                separation,... % Distance between nearest points of the wire (m)
                                vertSeparation) % Vertical separation of circuits (m). 0 for coplanar circuits
            function emfMat = calcEmfFromMatrix(msrCircuitAngle, drvCircuitAngle)
                % Deal with MatLab passing elements as arrays. Why? Why
                % make it so potentially convoluted for users?
                rows = size(msrCircuitAngle,1);
                cols = size(msrCircuitAngle,2);
                emfMat = zeros(rows, cols);
                for i = 1:rows
                    for j = 1:cols
                        rVec = InductiveFieldEmfForCircularCircuits.rVector(rDrvOuter, rMsrOuter, separation, vertSeparation, drvCircuitAngle(i,j), msrCircuitAngle(i,j));
                        accelVec = dIdt * [-sin(drvCircuitAngle(i,j)), cos(drvCircuitAngle(i,j)), 0];

                        emfMat(i,j) = rDrvOuter * rMsrOuter * dot(...
                            emfBetweenPoints(obj, accelVec, rVec),...
                            [-sin(msrCircuitAngle(i,j)), cos(msrCircuitAngle(i,j)) , 0]);
                    end
                end
            end
            
            emf = integral2(@(msrCircuitAngle, drvCircuitAngle) calcEmfFromMatrix(msrCircuitAngle, drvCircuitAngle), 0, 2 * pi, 0, 2 * pi);
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
                       