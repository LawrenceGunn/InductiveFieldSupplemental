function pass = RectangularCircuitTests()
    pass = true;
    pass = FaradayTests(pass);
    % pass = InductiveFieldTests(pass);
    if pass == true
        fprintf("All tests PASSED");
    else
        fprintf("Tests FAILED");
    end
end

function pass = InductiveFieldTests(currentPass)
    pass = currentPass;
    micro = 1000000.0;
    
    eField = InductiveFieldEmfForRectangularCircuits;

   % Test 13-15: wiresFromCorners test for single wire.
    wires21 = RectangularCircuitsCommon.wiresFromCorners({srcCorner1, srcCorner2}, srcOffset, false);
    EXPECT_NEAR(1, length(wires21), 0, "totalEmf test 14", pass);
    EXPECT_NEAR(1, wires21{1}.wireId, 0, "totalEmf test 14", pass);
    EXPECT_NEAR([0.0, 0.0507425, 0.0], wires21{1}.end, 0, "totalEmf test 14", pass);

    % Test 16-21: wiresFromCorners test for full driven wire loop.
    srcWires = RectangularCircuitsCommon.wiresFromCorners(srcCorners, srcOffset, true);
    EXPECT_NEAR(4, length(srcWires), 0, "totalEmf test 14", pass);
    EXPECT_NEAR(4, srcWires{4}.wireId, 0, "totalEmf test 14", pass);
    EXPECT_NEAR(srcCorner2, srcWires{1}.end, 0, "totalEmf test 14", pass);
    EXPECT_NEAR(srcCorner3, srcWires{3}.start, 0, "totalEmf test 14", pass);
    EXPECT_NEAR(srcCorner4, srcWires{4}.start, 0, "totalEmf test 14", pass);
    EXPECT_NEAR(srcCorner1, srcWires{4}.end, 0, "totalEmf test 14", pass);


end

function pass = FaradayTests(currentPass)
    pass = currentPass;
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
    pass = EXPECT_NEAR(-201.579, micro * emfFaradayAtClosest1, 0.001, "totalEmf test 14", pass);

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

