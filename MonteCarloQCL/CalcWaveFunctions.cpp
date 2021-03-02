#pragma warning(disable : 4996)

#include <math.h>
#include <complex> 
#include <fstream>
#include "Shoot.h"
#include "QCLMath.h"
#include "Constants.h"
#include "CalcWaveFunction.h"
#include "GenerateZSpace.h"

std::vector<WFStruct> CalculateWaveFunctions(std::vector<double> EigenEnergies, ZMaterialParmsStruct ZStruct)
{
    std::vector<WFStruct> WaveFunctions;

    //Temporary Struct for holding wavefunction during loop
    WFStruct TempWaveFunction;

    //Normalization Constant
    double Norm = 0;

    for (int k = 0; k < EigenEnergies.size(); k++)
    {
        //The spacing between grid points used to calculate the area under the |wavefunction|^2 for normalization  
        double DeltaZ;

        //Shoot wavefunctions using Shoot
        TempWaveFunction = Shoot(EigenEnergies[k], ZStruct);
        
        //Normalize Wavefunction
        for (int m = 0; m < TempWaveFunction.Wavefunction.size(); m++)
        {
            // Spacing between grid points
            DeltaZ = ZStruct.ZGridm[m + 1] - ZStruct.ZGridm[m];

            Norm += TempWaveFunction.Wavefunction[m] * TempWaveFunction.Wavefunction[m] * DeltaZ;
        }
        Norm = 1 / sqrt(Norm);

        TempWaveFunction.Wavefunction = TempWaveFunction.Wavefunction * Norm;

        //Add Normalized Wavefunction to vector of wavefunctions
        WaveFunctions.push_back(TempWaveFunction);

        Norm = 0;
    }
    
    return(WaveFunctions);
}

