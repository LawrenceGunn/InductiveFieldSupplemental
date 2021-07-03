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
        
        % This function calculates the EMF for one wire to another given the wire ends
        function emf = emfFromWireToWire(obj,...
                                         dIdt,...
                                         drvWireStart,...
                                         drvWireEnd,...
                                         msrWireStart,...
                                         msrWireEnd)
            % Calculate the slope m of the wires and the wire lengths
            drvLen = norm(drvWireEnd - drvWireStart);
            mDrv = (drvWireEnd - drvWireStart)/ drvLen;
            msrLen = norm(msrWireEnd - msrWireStart);
            mMsr = (msrWireEnd - msrWireStart)/ msrLen;
            dIdtVec = dIdt * mDrv;
            zeroVec = [0, 0, 0];
           
           
            function emfMat = calcEmfFromMatrix(sMsr, sDrv)
                % Deal with MatLab passing elements as arrays. Why? Why
                % make it so potentially convoluted for users?
                rows = size(sMsr,1);
                cols = size(sMsr,2);
                emfMat = zeros(rows, cols);
                for i = 1:rows
                    for j = 1:cols
                         rrVec = mMsr * sMsr(i,j) + msrWireStart - mDrv * sDrv(i,j) - drvWireStart;
                         emfMat(i,j) = dot(emfAtPoint(obj, dIdtVec, rrVec, zeroVec), mMsr);
                    end
                end
            end

            emf = integral2(@(sMsr, sDrv) calcEmfFromMatrix(sMsr, sDrv),...
                            0, msrLen, 0, drvLen,...
                            'AbsTol', 1e-12, 'RelTol', 1e-9);

        end
        
        % This function calculates the inline emf for two wire sets.
        function emf = inlineEmfForWires(obj,...
                                         dIdt,...
                                         drivenWire,...
                                         measuredWire)
            emf = emfFromWireToWire(obj,...
                                    dIdt,...
                                    drivenWire.start,...
                                    drivenWire.end,... 
                                    measuredWire.start,...
                                    measuredWire.end);
        end
        
        % This function calculates the inline emf for two wire sets.
        function inlineEmfs = inlineEmfForWireSetsOnWire(obj,...
                                                         dIdt,...
                                                         drivenWireSets,...
                                                         measuredWire)
        	emfTotal = 0.0;
            numDrivenWireSets = length(drivenWireSets);
            drivenWireEmfs = cell(numDrivenWireSets, 1);
            
            for i=1:numDrivenWireSets
                drivenWire = drivenWireSets{i};

                wireEmf = inlineEmfForWires(obj, dIdt, drivenWire, measuredWire);
                emfTotal = emfTotal + wireEmf;
                drivenWireEmfs{i}.drivenWireId = i;
                drivenWireEmfs{i}.emf = wireEmf;
            end
            
            inlineEmfs.emfTotal = emfTotal;
            inlineEmfs.drivenWireEmfs = drivenWireEmfs;
        end
        
        function obj = set.proposedAccelTermCoeff(obj, coeff)
            obj.proposedAccelTermCoeff = coeff;
        end
        
        function obj = set.conventionalAccelTermCoeff(obj, coeff)
            obj.conventionalAccelTermCoeff = coeff;
        end

        % This function calculates the inline emf for two sets of wires.
        function inlineEmfs = inlineEmfForWireSetsOnWireSets(obj,...
                                                             dIdt,...
                                                             drivenWireSets,...
                                                             measuredWireSets)
            emfTotalForSets = 0.0;
            numMeasuredWireSets = length(measuredWireSets);
            measuredWireEmfs = cell(numMeasuredWireSets, 1);

            for i=1:numMeasuredWireSets
                measuredWire = measuredWireSets{i};
                wireEmf = inlineEmfForWireSetsOnWire(obj, dIdt, drivenWireSets, measuredWire);
                emfTotalForSets = emfTotalForSets + wireEmf.emfTotal;
                measuredWireEmfs{i}.measuredWireId = measuredWire.wireId;
                measuredWireEmfs{i}.measuredWireEmfs = wireEmf;
            end

            inlineEmfs.emfTotal = emfTotalForSets;
            inlineEmfs.individualWireEmfs = measuredWireEmfs;
        end
    
    end
    
end
                       