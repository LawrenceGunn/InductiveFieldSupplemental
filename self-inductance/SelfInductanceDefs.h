/**********************************************************************

File     : StraightWireSelfInductance.h
Project  : Field Simulator
Purpose  : Header file for the finding numerical approximations of dtr/dt.
           2020/06/04

*/
#pragma once
#include <vector>
#include "Eigen/Dense"

namespace selfInductance {

    constexpr double PI = 3.1415926535897932384626433832795f;
    constexpr double VACUUM_PERMEABILITY = 1.2566370614e-6; // H/m (Henry's per meter) or N/A/A (Newtons per Ampere)
    constexpr double VACUUM_PERMITTIVITY = 8.85418781213e-12; // C^2 s^2/(kg m^2)
    constexpr double SPEED_OF_LIGHT = 299792458.0; // m/s
    constexpr double CHARGES_PER_COULOMB = 6241509074460762607.776;
    constexpr double HENRIES_TO_NANO_HENRIES = 1.0e9;
    constexpr double MICRO_HENRIES_TO_NANO_HENRIES = 1.0e3;
    constexpr double MILLIMETERS_TO_CENTIMETERS = 1.0e-1;
    constexpr double MILLIMETERS_TO_METERS = 1.0e-3;
    constexpr double UNINITIALIZED = std::numeric_limits<double>::quiet_NaN();

    using Vec = Eigen::Vector3d;
    using ForceVec = Vec;
    using Node = Eigen::Vector3d;
    using Polygon = Eigen::Vector3i;

    struct ChargeElement {
        Vec centroid;
        double area;
        double charge;
        double fractionOfCharge;
    };

    using VecList = std::vector<Vec>;
    using NodeList = std::vector<Node>;
    using PolygonList = std::vector<Polygon>;
    using ChargeElementList = std::vector<ChargeElement>;

    enum class InductanceTermsToUse {
        CONVENTIONAL_TERM,
        PROPOSED_TERM
    };
}
