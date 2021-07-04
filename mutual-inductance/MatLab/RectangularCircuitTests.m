function pass = RectangularCircuitTests(verbose)
    if ~exist('verbose','var')
        verbose = true;
    end
    
    pass = true;
    skipLongTests = false;
    
    testNumber = 1;
    FaradayTests();

    testNumber = 1;
    InductiveFieldTests();
    
    if pass == true
        fprintf("All tests PASSED");
    else
        fprintf("Tests FAILED");
    end

    function passed = InductiveFieldTests()
        fprintf("Running inductive field rectangular circuit tests\n");
        micro = 1000000.0;
        milli = 1000.0;

        eField = InductiveFieldEmfForRectangularCircuits;

        dIdtClosest = 6110.3534985431;

        srcCorner1 = [0.0, -0.0507425, 0.0];
        srcCorner2 = [0.0, 0.0507425, 0.0];
        srcCorner3 = [-0.101035, 0.0507425, 0.0];
        srcCorner4 = [-0.101035, -0.0507425, 0.0];
        srcCorners = {srcCorner1, srcCorner2, srcCorner3, srcCorner4};
        srcOffset = [0 0 0];

        msrCorner1 = [0.0, 0.0195775, 0.0];
        msrCorner2 = [0.0, -0.0195775, 0.0];
        msrCorner3 = [0.050515, -0.0195775, 0.0];
        msrCorner4 = [0.050515, 0.0195775, 0.0];
        msrCorners = {msrCorner1, msrCorner2, msrCorner3, msrCorner4};
        msrOffsetClosest = [0.000615, 0.0, 0.0];
        msrOffsetFarthest = [0.039175, 0.0, 0.0];

        % Test 1
        r1 = [2, 3, 0];
        a1 = [1, 1, 0];
        v1 = accelProduct(eField, a1, r1);
        unitR1 = r1/norm(r1);
        EXPECT_NEAR(micro * (cross(unitR1, cross(unitR1, a1)) - a1), micro * v1, 0.01, "accel product test 1");

        % Test 2
        v2 = accelProduct(eField, [1, 0, 0], [1, 1, 0]);
        EXPECT_NEAR(micro * [-1.5, 0.5, 0], micro * v2, 0.00001, "accel produc test 2");

        % Test 3
        v3 = accelProduct(eField, [1, 0, 0], [0, 1, 0]);
        EXPECT_NEAR(micro * [-2, 0, 0], micro * v3, 0.00001, "accel produc test 3");

        % Test 4: compare emf to old method. It seems like the old version
        % produced an EMF opposite of that calculated in the new code.
        r4Drv = [0.000615, 0.0, 0];
        r4Msr = [0.0, 0.0, 0];
        dIdt = [dIdtClosest, 0, 0]; % old used acceleration in x direction only*)
        eFieldV4 = emfAtPoint(eField, dIdt, r4Drv, r4Msr);
        EXPECT_NEAR([-0.993553, 0.0, 0.0],  eFieldV4, 0.000001, "emf at point test 4");

        % Tests 5-9: Get EMF at end points and middle point of Test 5
        % because the test is currently failing).
        dIdt5 = 1.0;
        dIdtVec5 = [dIdt5, 0.0, 0.0];
        srcStart5 = [-0.05, 0, 0];
        srcEnd5 = [0.05, 0, 0];
        msrPnt5 = [0.0, 0.002, 0];
        emf5 = emfAtPoint(eField, dIdtVec5, srcStart5, msrPnt5);
        EXPECT_NEAR([-2.00159, 0.0798084, 0.0], micro * emf5, 0.00001, "emf at point test 5");

        emf6 = emfAtPoint(eField, dIdtVec5, mean(srcStart5 + srcEnd5), msrPnt5);
        EXPECT_NEAR([-100, 0, 0], micro * emf6, 0.000001, "emf at point test 6");

        emf7 = emfAtPoint(eField, dIdtVec5, srcEnd5, msrPnt5);
        EXPECT_NEAR([-2.00159, -0.0798084, 0.0], micro * emf7, 0.00001, "emf at point test 7");
        EXPECT_NEAR(micro * emf5(1), micro * emf7(1), 0.000001, "emf at point test 8");
        EXPECT_NEAR(micro * emf5(2), -micro * emf7(2), 0.000001, "emf at point test 9");

        % Test 10: emfFromWireToWire
        srcStart11 = [0.0, -0.05, 0.0];
        srcEnd11 = [0.0, 0.05, 0.0];
        msrPnt11 = [0.002, 0.0, 0.0];
        msrStart13 = [-0.05, 0.002, 0.0];
        msrEnd13 = [0.05, 0.002, 0.0];
        emf10 = emfFromWireToWire(eField, dIdt5, srcStart5, srcEnd5, msrStart13, msrEnd13);
        EXPECT_NEAR(-0.0921054, micro * emf10, 0.000001, "wire to wire emf 10");

        % Test 11: emfFromWireToWire perpendicular to Test 13
        srcStart15 = [srcStart5(2), srcStart5(1), 0.0];
        srcEnd15 = [srcEnd5(2), srcEnd5(1), 0.0];
        msrStart15 = [msrStart13(2), msrStart13(1), 0.0];
        msrEnd15 = [msrEnd13(2), msrEnd13(1), 0.0];
        emf11 = emfFromWireToWire(eField, dIdt5, srcStart15, srcEnd15, msrStart15, msrEnd15);
        EXPECT_NEAR(-0.0921054, micro * emf11, 0.000001, "wire to wire emf 11");

        % Test 12: emfFromWireToWire as per Test 13 with measured wire at 30 degrees from parallel
        msrStart19 = [-0.025, 0.001, 0.0];
        msrEnd19 = [0.025, 0.001 + 0.05 * sin(deg2rad(30)), 0.0];
        emf12 = emfFromWireToWire(eField, dIdt5, srcStart5, srcEnd5, msrStart19, msrEnd19);
        EXPECT_NEAR(-0.0313341, micro * emf12, 0.000001, "wire to wire emf 12");

        % Test 13-15: wiresFromCorners test for single wire.
        wires21 = RectangularCircuitsCommon.wiresFromCorners({srcCorner1, srcCorner2}, srcOffset, false);
        EXPECT_NEAR(1, length(wires21), 0, "totalEmf test 13");
        EXPECT_NEAR(1, wires21{1}.wireId, 0, "totalEmf test 14");
        EXPECT_NEAR([0.0, 0.0507425, 0.0], wires21{1}.end, 0, "totalEmf test 15");

        % Test 16-21: wiresFromCorners test for full driven wire loop.
        srcWires = RectangularCircuitsCommon.wiresFromCorners(srcCorners, srcOffset, true);
        EXPECT_NEAR(4, length(srcWires), 0, "totalEmf test 16");
        EXPECT_NEAR(4, srcWires{4}.wireId, 0, "totalEmf test 17");
        EXPECT_NEAR(srcCorner2, srcWires{1}.end, 0, "totalEmf test 18");
        EXPECT_NEAR(srcCorner3, srcWires{3}.start, 0, "totalEmf test 19");
        EXPECT_NEAR(srcCorner4, srcWires{4}.start, 0, "totalEmf test 20");
        EXPECT_NEAR(srcCorner1, srcWires{4}.end, 0, "totalEmf test 21");

        % Test 22-27: wireFromCorners test with an offset.
        msrWiresClosest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetClosest, true);
        EXPECT_NEAR(4, length(msrWiresClosest), 0, "totalEmf test 22");
        EXPECT_NEAR(4, msrWiresClosest{4}.wireId, 0, "totalEmf test 23");
        EXPECT_NEAR(msrCorner2 + msrOffsetClosest, msrWiresClosest{1}.end, 0, "totalEmf test 24");
        EXPECT_NEAR(msrCorner3 + msrOffsetClosest, msrWiresClosest{3}.start, 0, "totalEmf test 25");
        EXPECT_NEAR(msrCorner4 + msrOffsetClosest, msrWiresClosest{4}.start, 0, "totalEmf test 26");
        EXPECT_NEAR(msrCorner1 + msrOffsetClosest, msrWiresClosest{4}.end, 0, "totalEmf test 27");

        % Test 28: inlineEmfForWires test for nearest wires at closest offset.
        if ~skipLongTests
            emf28 = inlineEmfForWires(eField, dIdtClosest, srcWires{1}, msrWiresClosest{1});
            EXPECT_NEAR(290.935, micro * emf28, 0.01, "inlineEmfForWires test 28");
        end

        % Test 29: inlineEmfForWires test for top driven and closest measured wires at closest offset.
        emf29 = inlineEmfForWires(eField, dIdtClosest, srcWires{2}, msrWiresClosest{1});
        EXPECT_NEAR(-13.3628, micro * emf29, 0.01, "inlineEmfForWires test 29");

        % Test 30: inlineEmfForWires test for top driven and top measured wires at closest offset.
        emf30 = inlineEmfForWires(eField, dIdtClosest, srcWires{2}, msrWiresClosest{4});
        EXPECT_NEAR(-55.2659, micro * emf30, 0.01, "inlineEmfForWires test 30");

        function emf = emfById(emfSet, id)
            for iDrv = 1:length(emfSet.drivenWireEmfs)
                if(emfSet.drivenWireEmfs{iDrv}.drivenWireId == id)
                    emf = emfSet.drivenWireEmfs{iDrv};
                    break;
                end
            end
        end

        % Tests 31-35: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf31 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{1});
        EXPECT_NEAR(290.929,  micro * emfById(emf31, 1).emf, 0.001, "inlineEmfForWireSetsOnWire test 31");
        EXPECT_NEAR(-13.3628, micro * emfById(emf31, 2).emf, 0.001, "inlineEmfForWireSetsOnWire test 32");
        EXPECT_NEAR(-44.007,  micro * emfById(emf31, 3).emf, 0.001, "inlineEmfForWireSetsOnWire test 33");
        EXPECT_NEAR(-13.3628, micro * emfById(emf31, 4).emf, 0.001, "inlineEmfForWireSetsOnWire test 34");
        if verbose
            fprintf("Starting wire to wire for old and proposed terms\n");
        end
        EXPECT_NEAR(220.196, micro * emf31.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 35");

        % Test 36: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf36 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{2});
        EXPECT_NEAR(-12.1687, micro * emf36.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 36");

        % Test 37: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf37 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{3});
        EXPECT_NEAR(-25.8672, micro * emf37.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 37");

        % Test 38: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf38 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{4});
        EXPECT_NEAR(-12.1687, micro * emf38.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 38");

        % Test 39: inlineEmfForWireSetsOnWireSets
        emf39 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires, msrWiresClosest);
        EXPECT_NEAR(169.992, micro * emf39.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 39");

        % Test 40-45: wireFromCorners test with the furthest offset.
        msrWiresFarthest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetFarthest, true);
        EXPECT_NEAR(4, length(msrWiresFarthest), 0, "wiresFromCorners test 40");
        EXPECT_NEAR(4, msrWiresFarthest{4}.wireId, 0, "wiresFromCorners test 41");
        EXPECT_NEAR(msrCorner2 + msrOffsetFarthest, msrWiresFarthest{1}.end, 0, "wiresFromCorners test 42");
        EXPECT_NEAR(msrCorner3 + msrOffsetFarthest, msrWiresFarthest{3}.start, 0, "wiresFromCorners test 43");
        EXPECT_NEAR(msrCorner4 + msrOffsetFarthest, msrWiresFarthest{4}.start, 0, "wiresFromCorners test 44");
        EXPECT_NEAR(msrCorner1 + msrOffsetFarthest, msrWiresFarthest{4}.end, 0, "wiresFromCorners test 45");

        % Test 46: inlineEmfForWireSetsOnWireSets
        emf46 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires, msrWiresFarthest);
        EXPECT_NEAR(12.3345, micro * emf46.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 46");

        experimentalXList = [0.000615, 0.000945, 0.002155, 0.003545, 0.005305,...
            0.006915, 0.009085, 0.011315, 0.014395, 0.017285, 0.020965,...
           0.026325, 0.031905, 0.039175];

        % Test 47: inlineEmfForWireOnWireSetsAtOffsets from experimental data
        if ~skipLongTests
            emfAtOffset47 = inlineEmfForWireOnWireSetsAtOffsets(eField, dIdtClosest, srcCorners, msrCorners, experimentalXList);
            if verbose
                fprintf("emfAtOffset47:\n");
                for i=1:length(emfAtOffset47)
                    fprintf("    x offset : %-10g   emf: %-10g\n", milli * emfAtOffset47{i}.xOffset, micro * emfAtOffset47{i}.emfTotal);
                end
            end
          EXPECT_NEAR(12.3345, micro * emfAtOffset47{14}.emfTotal, 0.0001, "inlineEmfForWireOnWireSetsAtOffsets test 47");
        end

        % Test 48: emfAtPoint for comparison to rectangular wires
        rDrv67 = [0, 0, 0];
        rMsr67 = [0.021, 0.02, 0];
        emf48 = emfAtPoint(eField, [0, dIdtClosest, 0], rDrv67, rMsr67);
        EXPECT_NEAR([0.0105226, -0.0321189, 0.0], emf48, 0.0001, "emfAtPoint test 48");

        % Test 49: inlineEmfForWireOnWireSetsAtOffsets from experimental data
        if ~skipLongTests
            emfAtOffset49 = inlineEmfForWireOnWireSetsAtOffsets(eField, dIdtClosest, srcCorners, msrCorners, experimentalXList);
            if verbose
                fprintf("emfAtOffset49:\n");
                for i=1:length(emfAtOffset49)
                    fprintf("    x offset : %-10g   emf: %-10g\n", milli * emfAtOffset49{i}.xOffset, micro * emfAtOffset49{i}.emfTotal);
                end
            end
            EXPECT_NEAR(12.3345, micro * emfAtOffset49{14}.emfTotal, 0.0001, "inlineEmfForWireOnWireSetsAtOffsets test 49");
        end

        % Run tests for old and proposed LW E-field terms individually.
        eField.conventionalAccelTermCoeff = 1;
        eField.proposedAccelTermCoeff = 0;

        % Tests 50-54
        emf50 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{1});
        if verbose
            fprintf("Front wire to all wires for existing double cross term only");
        end
        EXPECT_NEAR(47.8461, micro * emfById(emf50, 1).emf, 0.001, "inlineEmfForWireSetsOnWire test 50");
        EXPECT_NEAR(-13.3628,  micro * emfById(emf50, 2).emf, 0.001, "inlineEmfForWireSetsOnWire test 51");
        EXPECT_NEAR(-21.1205,  micro * emfById(emf50, 3).emf, 0.001, "inlineEmfForWireSetsOnWire test 52");
        EXPECT_NEAR(-13.3628,  micro * emfById(emf50, 4).emf, 0.001, "inlineEmfForWireSetsOnWire test 53");
        if verbose
            fprintf("Starting wire to wire for existing double cross term only");
        end
        EXPECT_NEAR(0, micro * emf50.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 54");
        
        % Test 55
        emf55 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{2});
        EXPECT_NEAR(0, micro * emf55.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 55");

        % Test 56
        emf56 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{3});
        EXPECT_NEAR(0, micro * emf56.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 56");

        % Test 57
        emf57 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{4});
        EXPECT_NEAR(0, micro * emf57.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 57");

        % Test 58
        emf58 = emf50.emfTotal + emf55.emfTotal + emf56.emfTotal + emf57.emfTotal;
        EXPECT_NEAR(0, micro * emf58, 0.001, "inlineEmfForWireSetsOnWire test 58");

        % Run tests for old and proposed LW E-field terms individually.
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;

        % Tests 59-63
        emf59 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{1});
        EXPECT_NEAR(243.083, micro * emfById(emf59, 1).emf, 0.001, "inlineEmfForWireSetsOnWire test 59");
        EXPECT_NEAR(0,  micro * emfById(emf59, 2).emf, 0.001, "inlineEmfForWireSetsOnWire test 60");
        EXPECT_NEAR(-22.8865, micro * emfById(emf59, 3).emf, 0.001, "inlineEmfForWireSetsOnWire test 61");
        EXPECT_NEAR(0, micro * emfById(emf59, 4).emf, 0.001, "inlineEmfForWireSetsOnWire test 62");
        if verbose
            fprintf("Starting wire to wire for proposed a term only");
        end
        EXPECT_NEAR(220.196, micro * emf59.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 63");

        % Test 64
        emf64 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{2});
        EXPECT_NEAR(-12.1687, micro * emf64.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 64");

        % Test 65
        emf65 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{3});
        EXPECT_NEAR(-25.8672, micro * emf65.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 65");

        % Test 66
        emf66 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{4});
        EXPECT_NEAR(-12.1687, micro * emf66.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 66");

        % Test 67
        emf67 = emf59.emfTotal + emf64.emfTotal + emf65.emfTotal +  emf66.emfTotal;
        EXPECT_NEAR(169.992, micro * emf67, 0.001, "inlineEmfForWireSetsOnWire test 67");
        
        % Test 68: wireset with vertical displacement
        srcWidth68 = 0.101035;
        msrWidth68 = 0.050515;
        sep68 = (srcWidth68 + msrWidth68)/2 + msrOffsetClosest(1);
        srcWires68 = RectangularCircuitsCommon.wiresFromDimensions([srcWidth68, 2 * 0.0507425], [0, 0, 0]);
        msrWires68 =RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [sep68, 0, 0]);

        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;

        emf68 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires68);
        EXPECT_NEAR(169.992, micro * emf68.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 68");

        % Test 69: wireset with vertical displacement
        msrWires69 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.0]);
        emf69 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires69);
        EXPECT_NEAR(-148.064, micro * emf69.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 69");

        % Test 70: wireset with vertical displacement
        msrWires70 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.01]);
        emf70 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires70);
        EXPECT_NEAR(-138.135, micro * emf70.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 70");

        % Test 71: wireset with vertical displacement
        msrWires71 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.02]);
        emf71 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires71);
        EXPECT_NEAR(-115.356, micro * emf71.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 71");

        % Test 72: wireset with vertical displacement
        msrWires72 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, -0.02]);
        emf72 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires72);
        EXPECT_NEAR(-115.356, micro * emf72.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 72");

        eField.proposedAccelTermCoeff = 0;
        eField.conventionalAccelTermCoeff = 1;

        % Test 73: wireset with vertical displacement
        msrWires73 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.0]);
        emf73 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires73);
        EXPECT_NEAR(0.0, micro * emf73.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 73");

        % Test 74: wireset with vertical displacement
        msrWires74 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.01]);
        emf74 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires74);
        EXPECT_NEAR(0.0, micro * emf74.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 74");

        % Test 75: wireset with vertical displacement
        msrWires75 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth68, 2 * 0.0195775], [0, 0, 0.02]);
        emf75 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires68, msrWires75);
        EXPECT_NEAR(0.0, micro * emf75.emfTotal, 0.001, "inlineEmfForWireSetsOnWireSets test 75");

        % Test 76: perpendicular wires
        eField.proposedAccelTermCoeff = 1;
        eField.conventionalAccelTermCoeff = 0;
        emf76 = emfFromWireToWire(eField, dIdtClosest, [1, 0, 0], [2, 0, 0], [2, 0.5, 0], [2, 2, 0]);
        EXPECT_NEAR(0.0, micro * emf76, 0.001, "emfFromWireToWire test 76");
        
        % Test 77: from comparison table non-zero perp wire issue
        sep77 = 0.005;
        srcWires77 = RectangularCircuitsCommon.wiresFromDimensions([0.1, 0.1], [0, 0, 0]);
        msrWires77 = RectangularCircuitsCommon.wiresFromDimensions([0.05, 0.04], [(0.1 + 0.05)/2 + sep77, 0, 0]);
        dot77 = dot(srcWires77{1}.end - srcWires77{1}.start, msrWires77{4}.start - msrWires77{4}.end);
        EXPECT_NEAR(0.0, dot77, 0.00001, "wiresFromDimensions test 77");

        % Test 78: from comparison table non-zero perp wire issue - emf being start/end - perp wires
        dIdt78 = 1000.0;
        emf78 = emfFromWireToWire(eField, dIdt78, srcWires77{1}.end, srcWires77{1}.start, msrWires77{4}.start, msrWires77{4}.end);
        EXPECT_NEAR(0.0, micro * emf78, 0.001, "emfFromWireToWire test 78");

        % Test 79: from comparison table non-zero perp wire issue - emf being start/end - nearest wires
        emf79 = emfFromWireToWire(eField, dIdt78, srcWires77{1}.end, srcWires77{1}.start, msrWires77{3}.end, msrWires77{3}.start);
        EXPECT_NEAR(23.7651, micro * emf79, 0.001, "emfFromWireToWire test 79");

        % Test 80: from comparison table non-zero perp wire issue - emf being start/end - nearest wires by index
        emf80 = inlineEmfForWires(eField, dIdt78, srcWires77{1}, msrWires77{3});
        EXPECT_NEAR(23.7651, micro * emf80, 0.001, "inlineEmfForWires test 80");

        % Test 81: from comparison table non-zero perp wire issue - emf being start/end - nearest wires by index
        emf81 = inlineEmfForWires(eField, dIdt78, srcWires77{1}, msrWires77{1});
        EXPECT_NEAR(-6.46048, micro * emf81, 0.001, "inlineEmfForWires test 81");

        % Test 82: from comparison table non-zero perp wire issue - emf being start/end - nearest wires by index
        emf82 = inlineEmfForWires(eField, dIdt78, srcWires77{1}, msrWires77{2});
        EXPECT_NEAR(0.0, micro * emf82, 0.001, "inlineEmfForWires test 82");
        
        passed = pass;
    end

    function passed = FaradayTests()
        fprintf("Running Faraday' law rectangular circuit tests\n");
        milli = 1000.0;
        micro = 1000000.0;
        emfFdy = FaradayEmfForRectangularCircuits;

        dIdtClosest = 6110.3534985431;

        srcCorner1 = [0.0, -0.0507425, 0.0];
        srcCorner2 = [0.0, 0.0507425, 0.0];
        srcCorner3 = [-0.101035, 0.0507425, 0.0];
        srcCorner4 = [-0.101035, -0.0507425, 0.0];
        srcCorners = {srcCorner1, srcCorner2, srcCorner3, srcCorner4};
        srcOffset = [0 0 0];

        msrCorner1 = [0.0, 0.0195775, 0.0];
        msrCorner2 = [0.0, -0.0195775, 0.0];
        msrCorner3 = [0.050515, -0.0195775, 0.0];
        msrCorner4 = [0.050515, 0.0195775, 0.0];
        msrCorners = {msrCorner1, msrCorner2, msrCorner3, msrCorner4};
        msrOffsetClosest = [0.000615, 0.0, 0.0];
        msrOffsetFarthest = [0.039175, 0.0, 0.0];

        srcWires = RectangularCircuitsCommon.wiresFromCorners(srcCorners, srcOffset, true);
        msrWiresClosest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetClosest, true);

        msrWiresFarthest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetFarthest, true);

        experimentalXList = [0.000615, 0.000945, 0.002155, 0.003545, 0.005305,...
            0.006915, 0.009085, 0.011315, 0.014395, 0.017285, 0.020965, ...
           0.026325, 0.031905, 0.039175];

        % Test 1: faradaysEmfWireToWireSet from experimental data for near
        emfFaradayAtClosest1 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{1}, msrWiresClosest);
        EXPECT_NEAR(-201.579, micro * emfFaradayAtClosest1, 0.001, "faradaysEmf test 1");

        % Test 2: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest2 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{2}, msrWiresClosest);
        EXPECT_NEAR(12.1687, micro * emfFaradayAtClosest2, 0.001, "faradaysEmf test 2");

        % Test 3: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest3 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{3}, msrWiresClosest);
        EXPECT_NEAR(7.24911, micro * emfFaradayAtClosest3, 0.001, "faradaysEmf test 3");

        % Test 4: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest4 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{4}, msrWiresClosest);
        EXPECT_NEAR(12.1687, micro * emfFaradayAtClosest4, 0.001, "faradaysEmf test 4");

        % Test 5: faradaysEmfWireToWireSet from experimental data for top near
        faradayClosestSum = emfFaradayAtClosest1 + emfFaradayAtClosest2 + emfFaradayAtClosest3 + emfFaradayAtClosest4;
        EXPECT_NEAR(-169.992, micro * faradayClosestSum, 0.001, "faradaysEmf test 5");

        % Test 6: faradaysEmfWireSetToWireSet from experimental data for closest separation
        faradayClosestFull = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires, msrWiresClosest);
        EXPECT_NEAR(-169.992, micro * faradayClosestFull.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 6");

        % Test 7: inlineEmfForWireOnWireSetsAtOffsets from experimental data
        if ~skipLongTests
            emfAtOffset7 = faradayEmfForWireOnWireSetsAtOffsets(emfFdy, dIdtClosest, srcCorners, msrCorners, experimentalXList);
            fprintf("emfsAtOffsets:\n");
            for i=1:length(emfAtOffset7)
                fprintf("    x offset : %-10g   emf: %-10g\n", milli * emfAtOffset7{i}.xOffset, micro * emfAtOffset7{i}.emf);
            end
            EXPECT_NEAR(-12.3345, micro * emfAtOffset7{14}.emf, 0.0001, "faradayEmfForWireOnWireSetsAtOffsets test 6");
        end

        % Test 8: wireset with vertical displacement
        srcWidth8 = 0.101035;
        msrWidth8 = 0.050515;
        sep8 = (srcWidth8 + msrWidth8)/2 + msrOffsetClosest(1);
        srcWires8 = RectangularCircuitsCommon.wiresFromDimensions([srcWidth8, (2 * 0.0507425)], [0, 0, 0]);
        msrWires8 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [sep8, 0, 0]);

        emf8 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires8);
        EXPECT_NEAR(-169.992, micro * emf8.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 8");

        % Test 9: wireset with vertical displacement
        msrWires9 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0]);
        emf9 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires9);
        EXPECT_NEAR(148.064, micro * emf9.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 9");

        % Test 10: wireset with vertical displacement
        msrWires10 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0.01]);
        emf10 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires10);
        EXPECT_NEAR(138.135, micro * emf10.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 10");

        % Test 11: wireset with vertical displacement
        msrWires11 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0.02]);
        emf11 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires11);
        EXPECT_NEAR(115.356, micro * emf11.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 11");

        % Test 12: wireset with vertical displacement
        msrWires12 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, -0.02]);
        emf12 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires12);
        EXPECT_NEAR(115.356, micro * emf12.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 12");

                
        passed = pass;
    end

    function passed = EXPECT_NEAR(target, actual, eps, msg)
        
        if(isvector(target) && (size(target,1) ~= 1 || size(target,2) ~= 1))
            if verbose
                fprintf("Test: %d   Expected: [%g,%g,%g]  Actual: [%g,%g,%g]\n", testNumber, target(1), target(2), target(3), actual(1), actual(2), actual(3));
            end

            diff = abs(target-actual);
            maxDiff = max(diff);
            if(maxDiff > eps)
                fprintf('Error for %s: Max difference is %f\n', msg, maxDiff);
                pass = false;
            end
        else
            if verbose
                fprintf("Test: %d   Expected: %g  Actual: %g\n", testNumber, target, actual);
            end
            if(isvector(target))
                targetScalar = target(1);
                actualScalar = actual(1);
            else
                targetScalar = target;
                actualScalar = actual;
            end
            diff = targetScalar-actualScalar;
            if(abs(diff) > eps)
                fprintf('Error for %s: Difference is %f\n', msg, diff);
                fprintf('  Target is: %f\n', targetScalar);
                fprintf('  Actual is: %f\n', actualScalar);
                pass = false;
            end
        end
        
        testNumber = testNumber + 1;
        passed = pass;
    end
end

