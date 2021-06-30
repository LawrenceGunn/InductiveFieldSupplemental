classdef (ConstructOnLoad) FaradayEmfForRectangularCircuits
    properties
        e0 = 8.85418782e-12
        c = 299792458
        m0
    end
    methods
        function obj = FaradayEmfForRectangularCircuits
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
           
        function emf = faradayEmf(obj,...
                                  dIdt,... % Change of current A/s
                                  rDrvOuter,... %Driven circuit radius (m)
                                  rMsrOuter,... % Measured circuit radius (m)
                                  separation,... % Distance between nearest points of the wire (m)
                                  vertSeparation) % Vertical separation of circuits (m). 0 for coplanar circuits
            function emfMat = calcEmfFromMatrix(rMsr, msrCircuitAngle, drvCircuitAngle)
                % Deal with MatLab passing elements as arrays. Why? Why
                % make it so potentially convoluted for users?
                rows = size(rMsr,1);
                cols = size(rMsr,2);
                emfMat = zeros(rows, cols);
                for i = 1:rows
                    for j = 1:cols
                        emfMat(i,j) = emfAtPoint(obj, dIdt, rMsr(i,j), rDrvOuter, rMsrOuter, separation, vertSeparation, msrCircuitAngle(i,j), drvCircuitAngle(i,j));
                    end
                end
            end
            
            emf = integral3(@(rMsr, msrCircuitAngle, drvCircuitAngle) calcEmfFromMatrix(rMsr, msrCircuitAngle, drvCircuitAngle), 0, rMsrOuter, 0, 2 * pi, 0, 2 * pi);
        end
           
    end
    methods (Static)
        function vec = rVector(rDrvOuter, ... % Driven circuit radius (m)
                               rMsrOuter,... % Measured circuit radius (m)
                               separation,... % Distance between nearest points of the wire (m)
                               vertSeparation,... % Vertical separation of circuits (m). 0 for coplanar circuits
                               rMsr,... % Measured circuit radius at evaluation point (m)
                               drvCircuitAngle,... % The angle around the driven circuit in radians
                               msrCircuitAngle) % The angle around the measured circuit in radians

                           vec = [
                               rDrvOuter * (1.0 - cos(drvCircuitAngle)) + separation + rMsrOuter + rMsr * cos(msrCircuitAngle),...
                               rMsr * sin(msrCircuitAngle) - rDrvOuter * sin(drvCircuitAngle),...
                               vertSeparation];
        end

    end
end
                       