/**********************************************************************

File     : MeshedStraightWireSelfInductance.h
Project  : Bach Field Simulator
Purpose  : Header file for the finding numerical approximations of dtr/dt.
           2020/06/04

*/
#pragma once

#include "SelfInductanceDefs.h"

namespace selfInductance {

    class MeshedStraightWireSelfInductance {
    public:
        MeshedStraightWireSelfInductance();

        /***
         * This sets up the internal variables and creates the base mesh template.
         * @param wireLength
         * @param wireRadius
         * @param numMeshRings The number of mesh rings to be used, typically 3 to 5.
         * @param layerThickness This is the thickness of a cylindrical layer of the wire. Note
         * that the wire length should be evenly divisible by the layer thickness, and if not
         * will be reduced to be so.
         * @param dIdt
         */
        void initialize(double wireLength,
                        double wireRadius,
                        int numMeshRings,
                        double layerThickness,
                        InductanceTermsToUse termsToUse,
                        double dIdt);

        std::tuple<ForceVec, ForceVec> forceBetweenElements(const ChargeElement& srcElement,
                                                            const ChargeElement& msrElement);

        std::tuple<ForceVec, ForceVec>
        forceBetweenLayers(const std::vector<ChargeElement>& srcElements,
                           const std::vector<ChargeElement>& msrElements);

        std::tuple<ForceVec, ForceVec> totalForceInWire();

        double inductanceInHenries();

        double inductanceInHenriesFromUnitForces();

        double layerThickness() const { return m_layerThickness; }

        const std::vector<ForceVec>& getForcesByLayerDistance();

        const std::vector<ForceVec>& getCumulativeForcesByLayer();

        const std::vector<ForceVec>& getSummedForceAtLayer();

        const ChargeElementList& getChargeMeshTemplate() const { return m_meshLayer0; }

        Vec accelerationVector() const { return m_accelerationVec; }

        int numberOfLayers() const { return m_numLayers; }

        int numNodesPerLayer() const { return m_numNodesPerLayer; }

    private:
        void createForceVsDistanceTable();

        void runCalcThread(int layerStart, int layerEnd);

        double m_wireLength{UNINITIALIZED};
        double m_wireRadius{UNINITIALIZED};
        double m_layerThickness{UNINITIALIZED};
        double m_dIdt{UNINITIALIZED};

        ChargeElementList m_meshLayer0;
        int m_numNodesPerLayer{0};
        int m_numLayers{0};
        double m_exactWireCrossSectionalArea{UNINITIALIZED};
        double m_wireCrossSectionalArea{UNINITIALIZED};
        double m_wireVolumeCubicMeters;
        double m_chargesPerCubicMeter{UNINITIALIZED};
        double m_coulombsPerCubicMeter{UNINITIALIZED};
        double m_totalChargeInWireInCoulombs{UNINITIALIZED};
        double m_chargePerLayerInCoulombs{UNINITIALIZED};
        Vec m_accelerationVec;
        Vec m_accelerationUnitVec;
        int m_coeffProposed{0};
        int m_coeffConventional{0};

        std::vector<ForceVec> m_forcesByLayerDistance;
        std::vector<ForceVec> m_cumulativeForcesByLayer;
        std::vector<ForceVec> m_summedForceAtLayer;
        std::vector<ForceVec> m_forcesByLayerDistanceUnit;
        std::vector<ForceVec> m_cumulativeForcesByLayerUnit;
        std::vector<ForceVec> m_summedForceAtLayerUnit;
        ForceVec m_totalForceInWire;
        ForceVec m_totalForceInWireUnit;

        double m_chargeDensity{UNINITIALIZED};
        double m_accelToForceCoeffPerChargeInCoulombsSqr{UNINITIALIZED};

        std::mutex m_calcMutex;
        std::condition_variable m_calcConditionVariable;
        std::atomic_int m_calcNumComplete;
        int m_calcNumThreads;
    };
}
