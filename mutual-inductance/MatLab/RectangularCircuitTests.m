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

    % Test 1
    dIdt1 = 6000;
    sep1 = 0.01;
    vertSep1 = 0.0;
    rDrv1 = 0.07;
    rMsr1 = 0.05;

    r1 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep1,  0, deg2rad(180));
    pass = EXPECT_NEAR([sep1, 0, 0], r1, 0.000001, "rVector test 1", pass);
    
    % Test 2
    r2 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep1,  0, 0);
    pass = EXPECT_NEAR([2 * rMsr1 + sep1, 0, 0], r2, 0.000001, "rVector test 2", pass);

    % Test 3
    r3 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep1, 0, deg2rad(90));
    pass = EXPECT_NEAR([ rMsr1 + sep1, rMsr1, 0], r3, 0.000001, "rVector test 3", pass);

    % Test 4
    r4 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep1, deg2rad(90), deg2rad(90));
    pass = EXPECT_NEAR([ rMsr1 + sep1 + rDrv1, rMsr1 - rDrv1, 0], r4, 0.000001, "rVector test 4", pass);

    % Test 5
    r5 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep1, deg2rad(180), 0);
    pass = EXPECT_NEAR([ 2 * rMsr1 + sep1 + 2 * rDrv1, 0, 0], r5, 0.000001, "rVector test 5", pass);

    % Test 6
    vertSep2 = 0.1;
    r6 = InductiveFieldEmfForRectangularCircuits.rVector(rDrv1, rMsr1, sep1, vertSep2, deg2rad(180), 0);
    pass = EXPECT_NEAR([ 2 * rMsr1 + sep1 + 2 * rDrv1, 0, vertSep2], r6, 0.000001, "rVector test 6", pass);

    % Test 7
    eField.conventionalAccelTermCoeff = 0;
    eField.proposedAccelTermCoeff = 1;
    emf7 = totalEmf(eField, dIdt1, rDrv1, rMsr1, sep1, vertSep1);
    pass = EXPECT_NEAR(76.4685, micro * emf7, 0.0001, "totalEmf test 7", pass);

    % Test 8 : EMF for dimensions matching experiments
    dIdt8 = 13000.0; % Change of current (A/s)
    radiusDrv8 = 0.05;% Radius of driven circuit (m)
    radiusMsr8 = 0.0375; % Radius of measured circuit (m)
    sep8 = 0.0075; % Separation (m) with experimental EMF approx 0.000125 V
    emf8 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, sep8, 0);
    pass = EXPECT_NEAR(120.294, micro * emf8, 0.001, "totalEmf test 8", pass);

    % Test 9 : EMF for dimensions matching other calculation
    dIdt9 = 6110.3534985431;
    rDrv9 = 0.05;
    rMsr9 = 0.02;
    oneMmSep9 = 0.001;
    emf9 = totalEmf(eField, dIdt9, rDrv9, rMsr9, oneMmSep9, 0);
    pass = EXPECT_NEAR(66.7189, micro * emf9, 0.002, "totalEmf test 9", pass);

    % Test 10 : EMF for dimensions matching experiments
    sep10 = 0.00075;
    emf10 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, sep10, 0);
    pass = EXPECT_NEAR(238.08, micro * emf10, 0.001, "totalEmf test 10", pass);

    % Test 11
    eField.conventionalAccelTermCoeff = 1;
    eField.proposedAccelTermCoeff = 0;
    emf11 = totalEmf(eField, dIdt1, rDrv1, rMsr1, sep1, vertSep1);
    pass = EXPECT_NEAR(0.0, micro * emf11, 0.0001, "totalEmf test 11", pass);

    % Test 12 : EMF for dimensions matching experiments
    emf12 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, sep8, 0);
    pass = EXPECT_NEAR(0.0, micro * emf12, 0.001, "totalEmf test 12", pass);

    % Test 13 : EMF for dimensions matching other calculation
    emf13 = totalEmf(eField, dIdt9, rDrv9, rMsr9, oneMmSep9, 0);
    pass = EXPECT_NEAR(0.0, micro * emf13, 0.002, "totalEmf test 13", pass);

    % Test 14 : EMF for dimensions matching experiments
    sep14 = 0.00075;
    emf14 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, sep14, 0);
    pass = EXPECT_NEAR(0.0, micro * emf14, 0.001, "totalEmf test 14", pass);

    % Test 15 : Inline  EMF for dimensions matching experiments
    vertSep15 = 0.00075;
    eField.conventionalAccelTermCoeff = 0;
    eField.proposedAccelTermCoeff = 1;
    emf15 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, -(radiusDrv8 + radiusMsr8), vertSep15);
    pass = EXPECT_NEAR(-966.781, micro * emf15, 0.001, "totalEmf test 15", pass);

    % Test 16 : EMF for dimensions matching experiments
    eField.conventionalAccelTermCoeff = 1;
    eField.proposedAccelTermCoeff = 0;
    emf16 = totalEmf(eField, dIdt8, radiusDrv8, radiusMsr8, -(radiusDrv8 + radiusMsr8), vertSep15);
    pass = EXPECT_NEAR(-0.0, micro * emf16, 0.001, "totalEmf test 16", pass);

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
    srcOffset = zeros(3,1);

    msrCorner1 = [0.0, 0.0195775, 0.0];
    msrCorner2 = [0.0, -0.0195775, 0.0];
    msrCorner3 = [0.050515, -0.0195775, 0.0];
    msrCorner4 = [0.050515, 0.0195775, 0.0];
    msrCorners = {msrCorner1, msrCorner2, msrCorner3, msrCorner4};
    msrOffsetClosest = [0.000615, 0.0, 0.0];
    msrOffsetFarthest = [0.039175, 0.0, 0.0];

    srcWires = wiresFromCorners(emfFdy, srcCorners, srcOffset, true);
    msrWiresClosest = wiresFromCorners(emfFdy, msrCorners, msrOffsetClosest, true);

    msrWiresFarthest = wiresFromCorners(emfFdy, msrCorners, msrOffsetFarthest, true);

    experimentalXList = [0.000615, 0.000945, 0.002155, 0.003545, 0.005305,...
        0.006915, 0.009085, 0.011315, 0.014395, 0.017285, 0.020965, ...
       0.026325, 0.031905, 0.039175];

    % Test 1: faradaysEmfWireToWireSet from experimental data for near
    emfFaradayAtClosest1 = faradaysEmfWireToWireSet(emfFdy, dIdtClosest, srcWires(1), msrWiresClosest);
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

