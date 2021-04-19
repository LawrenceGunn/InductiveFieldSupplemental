/**********************************************************************

File     : StraightWireSelfInductance.h
Project  : Field Simulator
Purpose  : Header file for the finding numerical approximations of dtr/dt.
           Revisions: Original definition by L. Gunn.
           2020/06/04

*/
#pragma once

#include "SelfInductanceDefs.h"
#include <cmath>

namespace selfInductance {

    class RosaStraightWireSelfInductance {
    public:
        static double inductance(double wireLengthMeters,
                                 double wireRadius,
                                 double muR = VACUUM_PERMEABILITY) {
            return VACUUM_PERMEABILITY / (2.0 * PI) * wireLengthMeters * (
                std::log((wireLengthMeters / wireRadius) * (1.0 + std::sqrt(1.0 + std::pow(wireRadius / wireLengthMeters, 2))))
                - std::sqrt(1.0 + std::pow(wireRadius / wireLengthMeters, 2))
                + muR / VACUUM_PERMEABILITY / 4.0
                + wireRadius / wireLengthMeters
            );
        }
    };
}
