function [results] = RectangularCircuitWireEmfs(separation,... % Separation in x in meters
                                                drvDimension,...
                                                msrDimension)
                            
    dIdt = 1000.0;
    micro = 1000000.0;
    results = writeWireToWireEmfsToTable(separation);
   
    function circuit = circuitAtSeparation(separation, proposed, conv)
        eField = InductiveFieldEmfForRectangularCircuits;
        eField.proposedAccelTermCoeff = proposed;
        eField.conventionalAccelTermCoeff = conv;

        srcWires = RectangularCircuitsCommon.wiresFromDimensions(drvDimension, [0,0,0]);
        msrOffsetX = (drvDimension(1) + msrDimension(1))/2;

        msrWires = RectangularCircuitsCommon.wiresFromDimensions(msrDimension, [msrOffsetX + separation, 0, 0]);
        circuit.eField = eField;
        circuit.src = srcWires;
        circuit.msr = msrWires;
    end

    function emfs = findMsrLoopEmfForSingleSrcWire(separation, srcWireId, proposed, conv)
        circuit = circuitAtSeparation(separation, proposed, conv);
        emfs = cell(4,1);

        for msrId=1:4
            emfs{msrId} = inlineEmfForWires(circuit.eField, dIdt, circuit.src{srcWireId}, circuit.msr{msrId});
        end
    end

    function emfs = findMsrWireEmfForSrcLoop(separation, msrId, proposed, conv)
        circuit = circuitAtSeparation(separation, proposed, conv);
        emfs = cell(4,1);

        for srcWireId=1:4
            emfs{srcWireId} = inlineEmfForWires(circuit.eField, dIdt, circuit.src{srcWireId}, circuit.msr{msrId});
        end
    end

    function result = writeWireToWireEmfsToTable(separation)
        fprintf("  Conventional accel term EMF's (uV)        Proposed term EMF's (uV)\n");
        fprintf("           i=1         2          3          4           Sum                    i=1        2           3         4          Sum\n");
   
        sumConvEmf = [0, 0, 0, 0];
        sumProposedEmf = [0, 0, 0, 0];
   
        for j = 1:4
            row = "j=1";
            if j ~= 1
                row = "  "+int2str(j);
            end
            
            convEmfs = findMsrWireEmfForSrcLoop(separation, j, 0, 1);
            convEmfRowTotal = 0.0;
            for i=1:4
                sumConvEmf(i) = sumConvEmf(i) + micro * convEmfs{i};
                convEmfRowTotal = convEmfRowTotal + micro * convEmfs{i};
            end
            
            leftCols = sprintf("%s     % -10g % -10g % -10g % -10g % -10g", row, micro * convEmfs{1}, micro * convEmfs{2}, micro * convEmfs{3}, micro * convEmfs{4}, convEmfRowTotal);

            proposedEmfs = findMsrWireEmfForSrcLoop(separation, j, 1, 0);            
            proposedEmfRowTotal = 0.0;
            for i=1:4
                sumProposedEmf(i) = sumProposedEmf(i) + micro * proposedEmfs{i};
                proposedEmfRowTotal = proposedEmfRowTotal + micro * proposedEmfs{i};
            end
            
            rightCols = sprintf("%s     % -10g % -10g % -10g % -10g % -10g", row, micro * proposedEmfs{1}, micro * proposedEmfs{2}, micro * proposedEmfs{3}, micro * proposedEmfs{4}, proposedEmfRowTotal);
            fprintf("%s    %s\n", leftCols, rightCols);

        end
        
        leftCols = sprintf("Sum     % -10g % -10g % -10g % -10g % -10g", sumConvEmf(1), sumConvEmf(2), sumConvEmf(3), sumConvEmf(4), sum(sumConvEmf));
        rightCols = sprintf("Sum    % -10g % -10g % -10g % -10g % -10g", sumProposedEmf(1), sumProposedEmf(2), sumProposedEmf(3), sumProposedEmf(4), sum(sumProposedEmf));
        fprintf("%s    %s\n", leftCols, rightCols);            

        result = true;
    end

end
