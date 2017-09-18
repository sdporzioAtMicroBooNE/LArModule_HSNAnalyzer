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

namespace AuxVertex
{

    // Decay vertex class and functions
    class DecayVertex
    {
    public:
        // Constructor and destructor
        DecayVertex();
        virtual ~DecayVertex();

        DecayVertex(double x, double y, double z, int parIdx1, int parIdx2, std::string parType1, std::string parType2, std::string direction1, std::string direction2);

        // Getters
        double GetX() const;
        double GetY() const;
        double GetZ() const;
        double GetParX(int par) const;
        double GetParY(int par) const;
        double GetParZ(int par) const;
        int GetParIdx1() const;
        int GetParIdx2() const;
        std::string GetParType1() const;
        std::string GetParType2() const;
        std::string GetDirection1() const;
        std::string GetDirection2() const;
        bool IsInsideTPC() const;
        bool IsDetLocAssigned() const;
        int GetChannelLoc(int plane) const;
        double GetTickLoc(int plane) const;
        int GetParChannelLoc(int par,int plane) const;
        double GetParTickLoc(int par,int plane) const;

        // Setters
        void SetChannelLoc(int channel0, int channel1, int channel2);
        void SetTickLoc(double tick0, double tick1, double tick2);
        void SetParChannelLoc(int par, int channel0, int channel1, int channel2);
        void SetParTickLoc(int par, double tick0, double tick1, double tick2);
        void SetParXYZ(int par, double x, double y, double z);
        void SetIsInsideTPC(bool val);
        void SetIsDetLocAssigned(bool val);

        // Printers
        void PrintInformation() const;

    private:
    bool fIsInsideTPC; // Whetehr the vertex is inside the TPC.
    bool fIsDetLocAssigned; // Whether channel/tick coordinates have been determined.
    double fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
    std::vector<double> fParX, fParY, fParZ; // Spatial coordinates of the parent of the vertex inside the detector.
    int fParIdx1,  fParIdx2; // Index of parent in track/shower vector (same for origin vertices).
    std::string fParType1, fParType2; // Type of parent ('t'=track,'s'=shower,'n'=neutral).
    std::string fDirection1, fDirection2; // Whether the vertex was coming from the origin ('start') or end ('end') of a track/shower.
    std::vector<int> fChannelLoc; // Nearest channel in each plane.
    std::vector<double> fTickLoc; // Nearest time tick in each plane.
    std::vector<std::vector<int>> fParChannelLoc; // Nearest channel in each plane for the vertex parent.
    std::vector<std::vector<double>> fParTickLoc; // Nearest time tick in each plane for the vertex parent.
};

double Distance(DecayVertex v1, DecayVertex v2);
DecayVertex MeanVertex(DecayVertex v1, DecayVertex v2);

} //END namespace AuxVertex

#endif