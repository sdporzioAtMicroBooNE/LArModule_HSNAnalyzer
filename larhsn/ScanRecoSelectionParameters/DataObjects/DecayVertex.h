/******************************************************************************
 * @file DecayVertex.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.cxx
 * ****************************************************************************/

#ifndef DECAYVERTEX_H
#define DECAYVERTEX_H

// C++ standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <exception>

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace AuxVertex
{

    // Decay vertex class and functions
    class DecayVertex
    {
    public:
        // Constructor and destructor
        DecayVertex();
        virtual ~DecayVertex();

        DecayVertex(double x, double y, double z, int parIdx1, int parIdx2, std::string parType1, std::string parType2);

        // Getters
        double GetX();
        double GetY();
        double GetZ();
        int GetParIdx1();
        int GetParIdx2();
        std::string GetParType1();
        std::string GetParType2();
        bool IsInsideTPC();
        int GetWireLoc(int plane);
        double GetTickLoc(int plane);

    private:
    bool fIsInsideTPC; // Is the vertex inside the TPC.
    double fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
    int fParIdx1,  fParIdx2; // Index of parent in track/shower vector (same for origin vertices).
    std::string fParType1, fParType2; // Type of parent ('t'=track,'s'=shower,'n'=neutral).
    std::vector<int> fWireLoc; // Nearest wire in each plane.
    std::vector<double> fTickLoc; // Nearest time tick in each plane.

    geo::GeometryCore const* fGeometry; // Pointer to the Geometry service
    detinfo::DetectorProperties const* fDetectorProperties; // Pointer to the Detector Properties
};

double Distance(DecayVertex v1, DecayVertex v2);
DecayVertex MeanVertex(DecayVertex v1, DecayVertex v2);


} //END namespace AuxVertex

#endif