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
        
        % Calculate the EMF from Faraday's law for one driven wire on the measured wire loop.
        function emf = faradaysEmfWireToWireSet(obj,...
                                                dIdt,...
                                                srcWire,...
                                                msrWireSet)
   srcWireStartFdy = srcWire["start"];
   srcWireVecFdy = srcWire["end"] - srcWire["start"];
   srcWireLengthFdy = Norm[srcWireVecFdy];
   srcWireNormalizedFdy = Normalize[srcWireVecFdy];
   dIdtVecFdy = dIdt srcWireNormalizedFdy;
   
   (* This algorithm assumes that the the measured wire set is a closed
   rectangle. The first two sides are used as perpendicular vectors. 
   Start
   by verifying this.  *)
   
   msrXWireStartFdy =  msrWireSet[[1]]["start"];
   msrXWireEndFdy =  msrWireSet[[1]]["end"];
   msrYWireStartFdy =  msrWireSet[[2]]["start"];
   msrYWireEndFdy =  msrWireSet[[2]]["end"];
   
   msrXWireVecFdy = msrXWireEndFdy - msrXWireStartFdy;
   msrXWireLengthFdy = Norm[msrXWireVecFdy];
   msrXWireNormalizedFdy = Normalize[msrXWireVecFdy];
   
   msrYWireVecFdy = msrYWireEndFdy - msrYWireStartFdy;
   msrYWireLengthFdy = Norm[msrYWireVecFdy];
   msrYWireNormalizedFdy = Normalize[msrYWireVecFdy];
   
   msrXYDotFdy = Dot[msrXWireVecFdy, msrYWireVecFdy]/msrXWireLengthFdy;
   If[Abs[msrXYDotFdy] > 0.0001, 
    Throw["Error in faradaysEmfWireToWireSet: Wires not \
perpendicular"]];
   
   mu0 /(4*Pi)
    NIntegrate[
     srcWirePosFdy = 
      srcWireStartFdy + srcWireLinearPosFdy srcWireNormalizedFdy;
     msrCircuitPosFdy = 
      msrXWireStartFdy + xmFdy msrXWireNormalizedFdy + 
       ymFdy msrYWireNormalizedFdy;
     Dot[
       Cross[dIdtVecFdy, msrCircuitPosFdy - srcWirePosFdy],
       {0, 0, 1}
       ]/Norm[msrCircuitPosFdy - srcWirePosFdy]^3,
      {srcWireLinearPosFdy, 0, srcWireLengthFdy},
     {xmFdy, 0,  msrXWireLengthFdy},
     {ymFdy, 0, msrYWireLengthFdy},
     Method -> "LocalAdaptive"]
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
                       