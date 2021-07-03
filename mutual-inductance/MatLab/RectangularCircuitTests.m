function pass = RectangularCircuitTests(verbose)
    if ~exist('verbose','var')
        verbose = true;
    end
    
    pass = true;
    skipLongTests = false;
    
    testNumber = 1;
    %pass = FaradayTests(pass, skipLongTests);

    testNumber = 1;
    pass = InductiveFieldTests(pass, skipLongTests);
    if pass == true
        fprintf("All tests PASSED");
    else
        fprintf("Tests FAILED");
    end

    function pass = InductiveFieldTests(currentPass, skipLongTests)
        pass = currentPass;
        micro = 1000000.0;

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
        pass = EXPECT_NEAR(micro * (cross(unitR1, cross(unitR1, a1)) - a1), micro * v1, 0.01, "accel product test 1", pass);

        % Test 2
        v2 = accelProduct(eField, [1, 0, 0], [1, 1, 0]);
        pass = EXPECT_NEAR(micro * [-1.5, 0.5, 0], micro * v2, 0.00001, "accel produc test 2", pass);

        % Test 3
        v3 = accelProduct(eField, [1, 0, 0], [0, 1, 0]);
        pass = EXPECT_NEAR(micro * [-2, 0, 0], micro * v3, 0.00001, "accel produc test 3", pass);

        % Test 4: compare emf to old method. It seems like the old version
        % produced an EMF opposite of that calculated in the new code.
        r4Drv = [0.000615, 0.0, 0];
        r4Msr = [0.0, 0.0, 0];
        dIdt = [dIdtClosest, 0, 0]; % old used acceleration in x direction only*)
        eFieldV4 = emfAtPoint(eField, dIdt, r4Drv, r4Msr);
        pass = EXPECT_NEAR([-0.993553, 0.0, 0.0],  eFieldV4, 0.000001, "emf at point test 4", pass);

        % Tests 5-9: Get EMF at end points and middle point of Test 5
        % because the test is currently failing).
        dIdt5 = 1.0;
        dIdtVec5 = [dIdt5, 0.0, 0.0];
        srcStart5 = [-0.05, 0, 0];
        srcEnd5 = [0.05, 0, 0];
        msrPnt5 = [0.0, 0.002, 0];
        emf5 = emfAtPoint(eField, dIdtVec5, srcStart5, msrPnt5);
        pass = EXPECT_NEAR([-2.00159, 0.0798084, 0.0], micro * emf5, 0.00001, "emf at point test 5", pass);

        emf6 = emfAtPoint(eField, dIdtVec5, mean(srcStart5 + srcEnd5), msrPnt5);
        pass = EXPECT_NEAR([-100, 0, 0], micro * emf6, 0.000001, "emf at point test 6", pass);

        emf7 = emfAtPoint(eField, dIdtVec5, srcEnd5, msrPnt5);
        pass = EXPECT_NEAR([-2.00159, -0.0798084, 0.0], micro * emf7, 0.00001, "emf at point test 7", pass);
        pass = EXPECT_NEAR(micro * emf5(1), micro * emf7(1), 0.000001, "emf at point test 8", pass);
        pass = EXPECT_NEAR(micro * emf5(2), -micro * emf7(2), 0.000001, "emf at point test 9", pass);

        % Test 10: emfFromWireToWire
        srcStart11 = [0.0, -0.05, 0.0];
        srcEnd11 = [0.0, 0.05, 0.0];
        msrPnt11 = [0.002, 0.0, 0.0];
        msrStart13 = [-0.05, 0.002, 0.0];
        msrEnd13 = [0.05, 0.002, 0.0];
        emf10 = emfFromWireToWire(eField, dIdt5, srcStart5, srcEnd5, msrStart13, msrEnd13);
        pass = EXPECT_NEAR(-0.0921054, micro * emf10, 0.000001, "wire to wire emf 10", pass);

        % Test 11: emfFromWireToWire perpendicular to Test 13
        srcStart15 = [srcStart5(2), srcStart5(1), 0.0];
        srcEnd15 = [srcEnd5(2), srcEnd5(1), 0.0];
        msrStart15 = [msrStart13(2), msrStart13(1), 0.0];
        msrEnd15 = [msrEnd13(2), msrEnd13(1), 0.0];
        emf11 = emfFromWireToWire(eField, dIdt5, srcStart15, srcEnd15, msrStart15, msrEnd15);
        pass = EXPECT_NEAR(-0.0921054, micro * emf11, 0.000001, "wire to wire emf 11", pass);

        % Test 12: emfFromWireToWire as per Test 13 with measured wire at 30 degrees from parallel
        msrStart19 = [-0.025, 0.001, 0.0];
        msrEnd19 = [0.025, 0.001 + 0.05 * sin(deg2rad(30)), 0.0];
        emf12 = emfFromWireToWire(eField, dIdt5, srcStart5, srcEnd5, msrStart19, msrEnd19);
        pass = EXPECT_NEAR(-0.0313341, micro * emf12, 0.000001, "wire to wire emf 12", pass);

        % Test 13-15: wiresFromCorners test for single wire.
        wires21 = RectangularCircuitsCommon.wiresFromCorners({srcCorner1, srcCorner2}, srcOffset, false);
        pass = EXPECT_NEAR(1, length(wires21), 0, "totalEmf test 13", pass);
        pass = EXPECT_NEAR(1, wires21{1}.wireId, 0, "totalEmf test 14", pass);
        pass = EXPECT_NEAR([0.0, 0.0507425, 0.0], wires21{1}.end, 0, "totalEmf test 15", pass);

        % Test 16-21: wiresFromCorners test for full driven wire loop.
        srcWires = RectangularCircuitsCommon.wiresFromCorners(srcCorners, srcOffset, true);
        pass = EXPECT_NEAR(4, length(srcWires), 0, "totalEmf test 16", pass);
        pass = EXPECT_NEAR(4, srcWires{4}.wireId, 0, "totalEmf test 17", pass);
        pass = EXPECT_NEAR(srcCorner2, srcWires{1}.end, 0, "totalEmf test 18", pass);
        pass = EXPECT_NEAR(srcCorner3, srcWires{3}.start, 0, "totalEmf test 19", pass);
        pass = EXPECT_NEAR(srcCorner4, srcWires{4}.start, 0, "totalEmf test 20", pass);
        pass = EXPECT_NEAR(srcCorner1, srcWires{4}.end, 0, "totalEmf test 21", pass);

        % Test 22-27: wireFromCorners test with an offset.
        msrWiresClosest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetClosest, true);
        pass = EXPECT_NEAR(4, length(msrWiresClosest), 0, "totalEmf test 22", pass);
        pass = EXPECT_NEAR(4, msrWiresClosest{4}.wireId, 0, "totalEmf test 23", pass);
        pass = EXPECT_NEAR(msrCorner2 + msrOffsetClosest, msrWiresClosest{1}.end, 0, "totalEmf test 24", pass);
        pass = EXPECT_NEAR(msrCorner3 + msrOffsetClosest, msrWiresClosest{3}.start, 0, "totalEmf test 25", pass);
        pass = EXPECT_NEAR(msrCorner4 + msrOffsetClosest, msrWiresClosest{4}.start, 0, "totalEmf test 26", pass);
        pass = EXPECT_NEAR(msrCorner1 + msrOffsetClosest, msrWiresClosest{4}.end, 0, "totalEmf test 27", pass);

        % Test 28: inlineEmfForWires test for nearest wires at closest offset.
        if ~skipLongTests
            emf28 = inlineEmfForWires(eField, dIdtClosest, srcWires{1}, msrWiresClosest{1});
            pass = EXPECT_NEAR(290.935, micro * emf28, 0.01, "inlineEmfForWires test 28", pass);
        end

        % Test 29: inlineEmfForWires test for top driven and closest measured wires at closest offset.
        emf29 = inlineEmfForWires(eField, dIdtClosest, srcWires{2}, msrWiresClosest{1});
        pass = EXPECT_NEAR(-13.3628, micro * emf29, 0.01, "inlineEmfForWires test 29", pass);

        % Test 30: inlineEmfForWires test for top driven and top measured wires at closest offset.
        emf30 = inlineEmfForWires(eField, dIdtClosest, srcWires{2}, msrWiresClosest{4});
        pass = EXPECT_NEAR(-55.2659, micro * emf30, 0.01, "inlineEmfForWires test 30", pass);

        function emf = emfById(emfSet, id)
            for i = 1:length(emfSet.drivenWireEmfs)
                if(emfSet.drivenWireEmfs{i}.drivenWireId == id)
                    emf = emfSet.drivenWireEmfs{i};
                    return;
                end
            end
        end

        % Tests 31-35: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf31 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{1});
        pass = EXPECT_NEAR(290.929,  micro * emfById(emf31, 1).emf, 0.001, "inlineEmfForWireSetsOnWire test 31", pass);
        pass = EXPECT_NEAR(-13.3628, micro * emfById(emf31, 2).emf, 0.001, "inlineEmfForWireSetsOnWire test 32", pass);
        pass = EXPECT_NEAR(-44.007,  micro * emfById(emf31, 3).emf, 0.001, "inlineEmfForWireSetsOnWire test 33", pass);
        pass = EXPECT_NEAR(-13.3628, micro * emfById(emf31, 4).emf, 0.001, "inlineEmfForWireSetsOnWire test 34", pass);
        if verbose
            fprintf("Starting wire to wire for old and proposed terms\n");
        end
        pass = EXPECT_NEAR(220.196, micro * emf31.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 35", pass);

        % Test 36: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf36 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{2});
        pass = EXPECT_NEAR(-12.1687, micro * emf36.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 36", pass);

        % Test 37: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf37 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{3});
        pass = EXPECT_NEAR(-25.8672, micro * emf37.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 37", pass);

        % Test 38: inlineEmfForWireSetsOnWire test for all driven wires and closest measured wire.
        emf38 = inlineEmfForWireSetsOnWire(eField, dIdtClosest, srcWires, msrWiresClosest{4});
        pass = EXPECT_NEAR(-12.1687, micro * emf38.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 38", pass);

        % Test 39: inlineEmfForWireSetsOnWireSets
        emf39 = inlineEmfForWireSetsOnWireSets(eField, dIdtClosest, srcWires, msrWiresClosest);
        pass = EXPECT_NEAR(169.992, micro * emf39.emfTotal, 0.001, "inlineEmfForWireSetsOnWire test 39", pass);

        % Test 40-45: wireFromCorners test with the furthest offset.
        msrWiresFarthest = RectangularCircuitsCommon.wiresFromCorners(msrCorners, msrOffsetFarthest, true);
        pass = EXPECT_NEAR(4, length(msrWiresFarthest), 0, "wiresFromCorners test 40", pass);
        pass = EXPECT_NEAR(4, msrWiresFarthest{4}.wireId, 0, "wiresFromCorners test 41", pass);
        pass = EXPECT_NEAR(msrCorner2 + msrOffsetFarthest, msrWiresFarthest{1}.end, 0, "wiresFromCorners test 42", pass);
        pass = EXPECT_NEAR(msrCorner3 + msrOffsetFarthest, msrWiresFarthest{3}.start, 0, "wiresFromCorners test 43", pass);
        pass = EXPECT_NEAR(msrCorner4 + msrOffsetFarthest, msrWiresFarthest{4}.start, 0, "wiresFromCorners test 44", pass);
        pass = EXPECT_NEAR(msrCorner1 + msrOffsetFarthest, msrWiresFarthest{4}.end, 0, "wiresFromCorners test 45", pass);

    end

    function pass = FaradayTests(currentPass, skipLongTests)
        pass = currentPass;
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
        pass = EXPECT_NEAR(-201.579, micro * emfFaradayAtClosest1, 0.001, "faradaysEmf test 1", pass);

        % Test 2: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest2 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{2}, msrWiresClosest);
        pass = EXPECT_NEAR(12.1687, micro * emfFaradayAtClosest2, 0.001, "faradaysEmf test 2", pass);

        % Test 3: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest3 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{3}, msrWiresClosest);
        pass = EXPECT_NEAR(7.24911, micro * emfFaradayAtClosest3, 0.001, "faradaysEmf test 3", pass);

        % Test 4: faradaysEmfWireToWireSet from experimental data for top near
        emfFaradayAtClosest4 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires{4}, msrWiresClosest);
        pass = EXPECT_NEAR(12.1687, micro * emfFaradayAtClosest4, 0.001, "faradaysEmf test 4", pass);

        % Test 5: faradaysEmfWireToWireSet from experimental data for top near
        faradayClosestSum = emfFaradayAtClosest1 + emfFaradayAtClosest2 + emfFaradayAtClosest3 + emfFaradayAtClosest4;
        pass = EXPECT_NEAR(-169.992, micro * faradayClosestSum, 0.001, "faradaysEmf test 5", pass);

        % Test 6: faradaysEmfWireSetToWireSet from experimental data for closest separation
        faradayClosestFull = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires, msrWiresClosest);
        pass = EXPECT_NEAR(-169.992, micro * faradayClosestFull.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 6", pass);

        % Test 7: inlineEmfForWireOnWireSetsAtOffsets from experimental data
        if ~skipLongTests
            emfAtOffset7 = faradayEmfForWireOnWireSetsAtOffsets(emfFdy, dIdtClosest, srcCorners, msrCorners, experimentalXList);
            fprintf("emfsAtOffsets:\n");
            for i=1:length(emfAtOffset7)
                fprintf("    x offset : %-10g   emf: %-10g\n", milli * emfAtOffset7{i}.xOffset, micro * emfAtOffset7{i}.emf);
            end
            pass = EXPECT_NEAR(-12.3345, micro * emfAtOffset7{14}.emf, 0.0001, "faradayEmfForWireOnWireSetsAtOffsets test 6", pass);
        end

        % Test 8: wireset with vertical displacement
        srcWidth8 = 0.101035;
        msrWidth8 = 0.050515;
        sep8 = (srcWidth8 + msrWidth8)/2 + msrOffsetClosest(1);
        srcWires8 = RectangularCircuitsCommon.wiresFromDimensions([srcWidth8, (2 * 0.0507425)], [0, 0, 0]);
        msrWires8 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [sep8, 0, 0]);

        emf8 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires8);
        pass = EXPECT_NEAR(-169.992, micro * emf8.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 8", pass);

        % Test 9: wireset with vertical displacement
        msrWires9 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0]);
        emf9 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires9);
        pass = EXPECT_NEAR(148.064, micro * emf9.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 9", pass);

        % Test 10: wireset with vertical displacement
        msrWires10 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0.01]);
        emf10 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires10);
        pass = EXPECT_NEAR(138.135, micro * emf10.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 10", pass);

        % Test 11: wireset with vertical displacement
        msrWires11 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, 0.02]);
        emf11 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires11);
        pass = EXPECT_NEAR(115.356, micro * emf11.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 11", pass);

        % Test 12: wireset with vertical displacement
        msrWires12 = RectangularCircuitsCommon.wiresFromDimensions([msrWidth8, (2 * 0.0195775)], [0, 0, -0.02]);
        emf12 = faradaysEmfWireSetToWireSet(emfFdy, dIdtClosest, srcWires8, msrWires12);
        pass = EXPECT_NEAR(115.356, micro * emf12.emfTotal, 0.001, "faradaysEmfWireSetToWireSet test 12", pass);

    end

    function passed = EXPECT_NEAR(target, actual, eps, msg, currentPass)
        passed = currentPass;
        
        if(isvector(target) && (size(target,1) ~= 1 || size(target,2) ~= 1))
            if verbose
                fprintf("Test: %d   Expected: [%g,%g,%g]  Actual: [%g,%g,%g]\n", testNumber, target(1), target(2), target(3), actual(1), actual(2), actual(3));
            end

            diff = abs(target-actual);
            maxDiff = max(diff);
            if(maxDiff > eps)
                fprintf('Error for %s: Max difference is %f\n', msg, maxDiff);
                passed = false;
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
                passed = false;
            end
        end
        testNumber = testNumber + 1;
    end
end

