#pragma once

#include "SelfInductanceDefs.h"

namespace selfInductance {

    class CircularMesh {
    public:

        /***
         * Constructor
         * @param radius In meters
         * @param numRings The number of mesh rings. Must be 1 or greater
         */
        CircularMesh(double radiusM, int numRings);

        const ChargeElementList& getCircularChargeElementList() const { return m_circularElementList; }

        const ChargeElementList& getQuadrantChargeElementList() const { return m_elementList; }

        const PolygonList getPolygonList() const { return m_polygonList; }

        const NodeList getNodeList() const { return m_nodeList; }

        double circularArea() const { return m_circularArea; }

    private:

        void generateCircularChargeElementList();

        void generateChargeElementList();

        void generatePolygonList();

        void generateNodeList();

        double m_radiusM{0.0};
        int m_numRings{0};

        ChargeElementList m_circularElementList;
        ChargeElementList m_elementList;
        PolygonList m_polygonList;
        NodeList m_nodeList;
        double m_circularArea{0.0};
    };
}
    