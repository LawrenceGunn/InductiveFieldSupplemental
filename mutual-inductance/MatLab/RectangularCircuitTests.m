function pass = RectangularCircuitTests()
    pass = true;
    skipLongTests = false;
    %pass = FaradayTests(pass, skipLongTests);
    pass = InductiveFieldTests(pass, skipLongTests);
    if pass == true
        fprintf("All tests PASSED");
    else
        fprintf("Tests FAILED");
    end
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
    emf5 = emfAtPoint(eField,dIdtVec5, srcStart5, msrPnt5);
    pass = EXPECT_NEAR([-2.00159, 0.0798084, 0.0], micro * emf5, 0.00001, "emf at point test 5", pass);

    emf6 = emfAtPoint(eField,dIdtVec5, mean(srcStart5 + srcEnd5), msrPnt5);
    pass = EXPECT_NEAR([-100, 0, 0], micro * emf6, 0.000001, "emf at point test 6", pass);

    emf7 = emfAtPoint(eField,dIdtVec5, srcEnd5, msrPnt5);
    pass = EXPECT_NEAR([-2.00159, -0.0798084, 0.0], micro * emf7, 0.00001, "emf at point test 7", pass);
    pass = EXPECT_NEAR(micro * emf5(1), micro * emf7(1), 0.000001, "emf at point test 8", pass);
    pass = EXPECT_NEAR(micro * emf5(2), -micro * emf7(2), 0.000001, "emf at point test 9", pass);
    
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
    if(isvector(target) && size(target,1) ~= 1 && size(target,2) ~= 1)
        diff = abs(target-actual);
        maxDiff = max(diff);
        if(maxDiff > eps)
            fprintf('Error for %s: Max difference is %f\n', msg, maxDiff);
            passed = false;
        end
    else
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
end

