#include "MeshedStraightWireSelfInductance.h"
#include "CircularMesh.h"
#include <thread>

using namespace selfInductance;

namespace {
    const double EST_CHARGE_VELOCITY = 0.001; // Assume a current velocity of 1 mm/s in the wire
    const double EST_CURRENT = 1.0; // Assume 1 amp
    const double COPPER_ATOMS_PER_CUBIC_METER = 8.4e28;
    const ForceVec ZER0_VEC = Vec(0, 0, 0);
}

MeshedStraightWireSelfInductance::MeshedStraightWireSelfInductance() {
}

void MeshedStraightWireSelfInductance::initialize(double wireLength,
                                                  double wireRadius,
                                                  int numMeshRings,
                                                  double layerThickness,
                                                  InductanceTermsToUse termsToUse,
                                                  double dIdt) {
    m_wireLength = wireLength;
    m_wireRadius = wireRadius;
    m_numLayers = static_cast<int>(std::round(wireLength / layerThickness));
    m_layerThickness = wireLength / m_numLayers;

    m_dIdt = dIdt;

    m_coeffProposed = termsToUse == InductanceTermsToUse::CONVENTIONAL_TERM ? 0 : 1;
    m_coeffConventional = termsToUse == InductanceTermsToUse::PROPOSED_TERM ? 0 : 1;

    CircularMesh meshBuilder(m_wireRadius, numMeshRings);
    m_meshLayer0 = meshBuilder.getCircularChargeElementList();
    m_numNodesPerLayer = m_meshLayer0.size();

    m_exactWireCrossSectionalArea = PI * wireRadius * wireRadius;
    m_wireCrossSectionalArea = meshBuilder.circularArea();

    m_wireVolumeCubicMeters = m_exactWireCrossSectionalArea * m_wireLength;

    m_chargesPerCubicMeter = COPPER_ATOMS_PER_CUBIC_METER;
    m_coulombsPerCubicMeter = m_chargesPerCubicMeter / CHARGES_PER_COULOMB;

    m_totalChargeInWireInCoulombs = m_coulombsPerCubicMeter * m_wireVolumeCubicMeters;

    // At this point the total charge in the wire is known. Now it must be divided
    // up into each charge volume. To start this, determine the charge per layer.
    m_chargePerLayerInCoulombs = m_totalChargeInWireInCoulombs / m_numLayers;

    // Now take that charge and distribute it within the layer mesh template. Each
    // element takes a fraction of its charge which is defined by its area divided by
    // the total mesh area.
    for (auto& element : m_meshLayer0) {
        double fraction = element.area / m_wireCrossSectionalArea;
        double elementChargeInCoulombs = m_chargePerLayerInCoulombs * fraction;
        element.charge = elementChargeInCoulombs;
        element.fractionOfCharge = fraction / m_numLayers;
        element.centroid.z() = 0.5 * m_layerThickness;
    }

    double accelZ = m_dIdt / (meshBuilder.circularArea() * m_coulombsPerCubicMeter);

    m_accelerationVec = {0, 0, accelZ};
    m_accelerationUnitVec = {0, 0, 1.0};

    double c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    m_accelToForceCoeffPerChargeInCoulombsSqr = 1.0 / (4.0 * PI * c2 * VACUUM_PERMITTIVITY);
}

std::tuple<ForceVec, ForceVec>
MeshedStraightWireSelfInductance::forceBetweenElements(const ChargeElement& srcElement,
                                                       const ChargeElement& msrElement) {
    Vec rVec = msrElement.centroid - srcElement.centroid;
    double r = rVec.norm();
    if (std::abs(r) > std::numeric_limits<double>::epsilon()) {
        double rInv = 1.0 / r;
        ForceVec accelComponents = -m_coeffProposed * (m_accelerationVec * rInv);
        ForceVec accelComponentsUnit = -m_coeffProposed * (m_accelerationUnitVec * rInv);
        if (m_coeffConventional != 0) {
            accelComponents += rVec.cross(rVec.cross(m_accelerationVec)) * (rInv * rInv * rInv);
            accelComponentsUnit += rVec.cross(rVec.cross(m_accelerationUnitVec)) * (rInv * rInv * rInv);
        }
        return std::make_tuple(
            m_accelToForceCoeffPerChargeInCoulombsSqr * srcElement.charge * msrElement.charge * accelComponents,
            srcElement.fractionOfCharge * msrElement.fractionOfCharge * accelComponentsUnit);
    } else {
        return std::make_tuple(
            ZER0_VEC,
            ZER0_VEC);
    }
}

std::tuple<ForceVec, ForceVec>
MeshedStraightWireSelfInductance::forceBetweenLayers(const std::vector<ChargeElement>& srcElements,
                                                     const std::vector<ChargeElement>& msrElements) {
    ForceVec forceSum(0.0, 0.0, 0.0);
    ForceVec forceUnitSum(0.0, 0.0, 0.0);
    for (int iMsr = 0; iMsr < msrElements.size(); ++iMsr) {
        for (int iSrc = 0; iSrc < srcElements.size(); ++iSrc) {
            auto[force, forceUnit] = forceBetweenElements(srcElements[iSrc], msrElements[iMsr]);
            forceSum += force;
            forceUnitSum += forceUnit;
        }
    }

    return std::make_tuple(forceSum, forceUnitSum);
}

std::tuple<ForceVec, ForceVec>
MeshedStraightWireSelfInductance::totalForceInWire() {

    if (m_forcesByLayerDistance.empty()) {
        createForceVsDistanceTable();
    }

    m_totalForceInWire.setConstant(0.0);
    m_totalForceInWireUnit.setConstant(0.0);

    m_summedForceAtLayer.resize(m_numLayers);
    m_summedForceAtLayerUnit.resize(m_numLayers);

    for (int iLayer = 0; iLayer < m_numLayers; ++iLayer) {
        // The force for a layer will always include the self inductance for that layer.
        ForceVec forceForLayer = m_cumulativeForcesByLayer[0];
        ForceVec forceForLayerUnit = m_cumulativeForcesByLayerUnit[0];

        // For the low index side of the wire...
        int numNodesLow = iLayer;
        if (numNodesLow > 0) {
            forceForLayer += m_cumulativeForcesByLayer[numNodesLow];
            forceForLayerUnit += m_cumulativeForcesByLayerUnit[numNodesLow];
        }

        int numNodesHigh = m_numLayers - iLayer - 1;
        if (numNodesHigh > 0) {
            forceForLayer += m_cumulativeForcesByLayer[numNodesHigh];
            forceForLayerUnit += m_cumulativeForcesByLayerUnit[numNodesHigh];
        }

        m_summedForceAtLayer[iLayer] = forceForLayer;
        m_summedForceAtLayerUnit[iLayer] = forceForLayerUnit;

        m_totalForceInWire += forceForLayer;
        m_totalForceInWireUnit += forceForLayerUnit;
    }

    return std::make_tuple(m_totalForceInWire, m_totalForceInWireUnit);
}

double MeshedStraightWireSelfInductance::inductanceInHenries() {
    if (m_summedForceAtLayer.empty()) {
        totalForceInWire();
    }

    double areaSqr = m_exactWireCrossSectionalArea * m_exactWireCrossSectionalArea;
    double chargeDensityInCoulombsPerCubicMeter2 = m_coulombsPerCubicMeter * m_coulombsPerCubicMeter;
    return m_totalForceInWire(2) / (chargeDensityInCoulombsPerCubicMeter2 * areaSqr * m_accelerationVec(2));
}

double MeshedStraightWireSelfInductance::inductanceInHenriesFromUnitForces() {
    if (m_summedForceAtLayerUnit.empty()) {
        totalForceInWire();
    }

    double forceZ = m_totalForceInWireUnit(2);
    double basicCoeff = 1.0/(4.0 * PI * VACUUM_PERMITTIVITY * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    double lengthSqr = m_wireLength * m_wireLength;
    double inductance = basicCoeff * forceZ * lengthSqr;

    //double areaSqr = m_exactWireCrossSectionalArea * m_exactWireCrossSectionalArea;
    //double chargeDensityInCoulombsPerCubicMeter2 = m_coulombsPerCubicMeter * m_coulombsPerCubicMeter;
    return inductance;// / (chargeDensityInCoulombsPerCubicMeter2 * areaSqr * m_accelerationVec(2));
}

void MeshedStraightWireSelfInductance::createForceVsDistanceTable() {
//    VecList layerZero = generateNodesForLayer(0);
//    VecList layerN = generateNodesForLayer(0);
    m_forcesByLayerDistance.resize(m_numLayers);
    m_forcesByLayerDistanceUnit.resize(m_numLayers);

    m_calcNumComplete = 0;

    m_calcNumThreads = std::thread::hardware_concurrency();

    int numLayersPerThread = m_numLayers / (double) m_calcNumThreads;
    int nextStartLayer = 0;
    std::vector<std::thread> threads(m_calcNumThreads);
    for (int iThread = 0; iThread < m_calcNumThreads; ++iThread) {
        int startLayer = nextStartLayer;
        int endLayer = startLayer + numLayersPerThread - 1;
        if (iThread == m_calcNumThreads - 1) {
            endLayer = m_numLayers - 1;
        }

        threads[iThread] = std::thread([this, startLayer, endLayer] { runCalcThread(startLayer, endLayer); });
        nextStartLayer = endLayer + 1;
    }

    std::unique_lock<std::mutex> lock(m_calcMutex);
    m_calcConditionVariable.wait(lock, [this] {
        return m_calcNumThreads == m_calcNumComplete;
    });

    for (auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    m_cumulativeForcesByLayer.resize(m_numLayers);
    m_cumulativeForcesByLayerUnit.resize(m_numLayers);

    m_cumulativeForcesByLayer[0] = m_forcesByLayerDistance[0];
    m_cumulativeForcesByLayerUnit[0] = m_forcesByLayerDistanceUnit[0];

    if (m_numLayers > 1) {
        m_cumulativeForcesByLayer[1] = m_forcesByLayerDistance[1];
        m_cumulativeForcesByLayerUnit[1] = m_forcesByLayerDistanceUnit[1];

        for (int iLayer = 2; iLayer < m_numLayers; ++iLayer) {
            m_cumulativeForcesByLayer[iLayer] =
                m_cumulativeForcesByLayer[iLayer - 1] +
                m_forcesByLayerDistance[iLayer];

            m_cumulativeForcesByLayerUnit[iLayer] =
                m_cumulativeForcesByLayerUnit[iLayer - 1] +
                m_forcesByLayerDistanceUnit[iLayer];
        }
    }
}

void MeshedStraightWireSelfInductance::runCalcThread(int layerStart, int layerEnd) {
    ChargeElementList meshLayerN = m_meshLayer0;

    double centroidHeight0 = 0.5 * m_layerThickness;
    for (int iLayer = layerStart; iLayer <= layerEnd; ++iLayer) {

        // Increment the k (aka z) value of the node
        for (auto& element : meshLayerN) {
            element.centroid.z() = centroidHeight0 + iLayer * m_layerThickness;
        }

        auto[force, forceUnit] = forceBetweenLayers(meshLayerN, m_meshLayer0);
        m_forcesByLayerDistance[iLayer] = force;
        m_forcesByLayerDistanceUnit[iLayer] = forceUnit;
    }

    int currentCount = ++m_calcNumComplete;

    if (m_calcNumComplete == m_calcNumThreads) {
        m_calcConditionVariable.notify_all();
    }
}

const std::vector<ForceVec>& MeshedStraightWireSelfInductance::getForcesByLayerDistance() {
    if (m_forcesByLayerDistance.empty()) {
        createForceVsDistanceTable();
    }
    return m_forcesByLayerDistance;
}

const std::vector<ForceVec>& MeshedStraightWireSelfInductance::getCumulativeForcesByLayer() {
    if (m_cumulativeForcesByLayer.empty()) {
        createForceVsDistanceTable();
    }
    return m_cumulativeForcesByLayer;
}

const std::vector<ForceVec>& MeshedStraightWireSelfInductance::getSummedForceAtLayer() {
    return m_summedForceAtLayer;
}
