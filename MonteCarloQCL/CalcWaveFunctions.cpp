#pragma warning(disable : 4996)

#include <math.h>
#include <complex> 
#include <fstream>
#include "Shoot.h"
#include "QCLMath.h"
#include "Constants.h"
#include "CalcWaveFunction.h"

std::vector<WFStruct> CalculateWaveFunctions(std::vector<double> EigenEnergies, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms)
{
    std::vector<WFStruct> WaveFunctions;

    for (int k = 0; k < EigenEnergies.size(); k++)
    {
        WaveFunctions.push_back(Shoot(EigenEnergies[k], ZGridm, Potential, ms));
    }

    return(WaveFunctions);
}

