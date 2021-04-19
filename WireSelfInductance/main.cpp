#include "RosaStraightWireSelfInductance.h"
#include "MeshedStraightWireSelfInductance.h"
#include <iostream>
#include "args.hxx"

using namespace selfInductance;

int main(int argc,
         char** argv) {

    args::ArgumentParser parser("Self-inductance solver for a straight wire.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<double> length(parser,
                                   "length-mm",
                                   "Length of the wire in millimeters",
                                   {'l', "wire-length-mm"});
    args::ValueFlag<double> diameter(parser,
                                     "diameter-mm",
                                     "Diameter of the wire in millimeters",
                                     {'d', "wire-diameter-mm"});
    args::ValueFlag<int> numMeshRings(parser,
                                      "mesh-rings",
                                      "Number of mesh element rings",
                                      {'n', "num-mesh-rings"});
    args::ValueFlag<double> layerThicknessArg(parser,
                                              "layer-thickness",
                                              "Thickness of mesh elements",
                                              {'t', "layer-thickness-mm"});

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cout << parser << std::endl;
        return 0;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl << parser << std::endl;
        return 1;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl << parser << std::endl;
        return 1;
    }

    try {
        double wireLengthMm = args::get(length);
        double wireDiameterMm = args::get(diameter);
        int meshRings = args::get(numMeshRings);
        double layerThicknessMm = args::get(layerThicknessArg);

        if (meshRings < 1) {
            std::cerr << "ERROR: The number of mesh rings must be greater than 0" << std::endl;
            return 1;
        }

        double rosaInductanceMicroHenries = RosaStraightWireSelfInductance::inductance(
            MILLIMETERS_TO_METERS * wireLengthMm,
            MILLIMETERS_TO_METERS * 0.5 * wireDiameterMm);

        MeshedStraightWireSelfInductance proposedSolver;
        proposedSolver.initialize(MILLIMETERS_TO_METERS * wireLengthMm,
                                       MILLIMETERS_TO_METERS * 0.5 * wireDiameterMm,
                                       meshRings,
                                       MILLIMETERS_TO_METERS * layerThicknessMm,
                                       InductanceTermsToUse::PROPOSED_TERM,
                                       1.0);
        double proposedInductanceInHenries = proposedSolver.inductanceInHenries();

        MeshedStraightWireSelfInductance conventionalSolver;
        conventionalSolver.initialize(MILLIMETERS_TO_METERS * wireLengthMm,
                                       MILLIMETERS_TO_METERS * 0.5 * wireDiameterMm,
                                       meshRings,
                                       MILLIMETERS_TO_METERS * layerThicknessMm,
                                       InductanceTermsToUse::CONVENTIONAL_TERM,
                                       1.0);
        double conventionalInductanceInHenries = conventionalSolver.inductanceInHenries();

        std::cout << "Rosa inductance (nH): " << HENRIES_TO_NANO_HENRIES * rosaInductanceMicroHenries
                  << std::endl;
        std::cout << "Proposed inductance (nH): " << -HENRIES_TO_NANO_HENRIES * proposedInductanceInHenries
                  << std::endl;
        std::cout << "Conventional inductance (nH): " << -HENRIES_TO_NANO_HENRIES * conventionalInductanceInHenries
                  << std::endl;
    }
    catch (const std::exception& ex) {
        std::cerr << "EXCEPTION: " << ex.what() << std::endl;
    }

    return 0;
}
