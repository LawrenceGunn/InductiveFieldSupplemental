//
// Created on 7/4/20.
//

#include "CircularMesh.h"

using namespace selfInductance;

namespace {
    constexpr double INTERSECTION_TIME_RECALC_EPS = 1.0e-6;
    constexpr double ZERO_HEIGHT = 0.0;
    constexpr double ZERO_CHARGE = 0.0;
}

CircularMesh::CircularMesh(double radiusM, int numRings) :
    m_radiusM{radiusM},
    m_numRings{numRings} {
    generateCircularChargeElementList();
}

void CircularMesh::generateCircularChargeElementList() {
    generateChargeElementList();

    m_circularArea = 0.0;

    for(const auto& element : m_elementList) {
        m_circularElementList.push_back(element);
    }

    for (const auto& charge : m_elementList) {
        m_circularArea += 4.0*charge.area;

        ChargeElement charge2 = charge;
        charge2.centroid.x() = -charge2.centroid.x();
        m_circularElementList.push_back(charge2);

        ChargeElement charge3 = charge;
        charge3.centroid.x() = -charge3.centroid.x();
        charge3.centroid.y() = -charge3.centroid.y();
        m_circularElementList.push_back(charge3);

        ChargeElement charge4 = charge;
        charge4.centroid.y() = -charge4.centroid.y();
        m_circularElementList.push_back(charge4);
    }
}

void CircularMesh::generateChargeElementList() {
    generatePolygonList();

    for (const auto& polygon : m_polygonList) {
        double xa = m_nodeList[polygon(0)].x();
        double ya = m_nodeList[polygon(0)].y();
        double xb = m_nodeList[polygon(1)].x();
        double yb = m_nodeList[polygon(1)].y();
        double xc = m_nodeList[polygon(2)].x();
        double yc = m_nodeList[polygon(2)].y();
        Vec centroid((xa+xb+xc)/3.0, (ya+yb+yc)/3.0, ZERO_HEIGHT);
        double area = std::abs((xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb))/2.0);
        m_elementList.push_back({centroid, area, ZERO_CHARGE});
    }
}

void CircularMesh::generatePolygonList() {
    generateNodeList();

    int iLow = 1;
    int iHigh = 2;
    int numPolygonsForRing = 1;
    for (int iRing = 1; iRing <= m_numRings; ++iRing) {

        int iRingLow = iLow;
        int iRingHigh = iHigh;
        int numPolyAdded = 0;
        while (true) {
            m_polygonList.push_back({iRingLow - 1, iRingHigh - 1, iRingHigh});
            ++numPolyAdded;

            if (numPolyAdded >= numPolygonsForRing) {
                break;
            }

            ++iRingHigh;
            m_polygonList.push_back({iRingLow - 1, iRingHigh - 1, iRingLow});
            ++numPolyAdded;

            if (numPolyAdded >= numPolygonsForRing) {
                break;
            }

            ++iRingLow;
        }
        iLow = iHigh;
        iHigh += iRing + 1;
        numPolygonsForRing += 2;
    }
}

void CircularMesh::generateNodeList() {
    m_nodeList.push_back({0.0, 0.0, 0.0});
    for (int iRing = 1; iRing <= m_numRings; ++iRing) {
        double ringRadius = m_radiusM * iRing / static_cast<double>(m_numRings);
        for (int iNode = 0; iNode <= iRing; ++iNode) {
            m_nodeList.push_back({ringRadius * std::cos(0.5 * PI * iNode / static_cast<double>(iRing)),
                                  ringRadius * std::sin(0.5 * PI * iNode / static_cast<double>(iRing)),
                                  0.0}
            );
        }
    }
}
