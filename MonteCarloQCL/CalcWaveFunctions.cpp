#pragma warning(disable : 4996)

#include <math.h>
#include <complex> 
#include <fstream>
#include "Shoot.h"
#include "QCLMath.h"
#include "Constants.h"
#include "CalcWaveFunction.h"
#include "GenerateZSpace.h"

std::vector<WFStruct> CalculateWaveFunctions(std::vector<double> EigenEnergies, ZMaterialParmsStruct ZStruct, std::vector<double> Potential)
{
    std::vector<WFStruct> WaveFunctions;

    for (int k = 0; k < EigenEnergies.size(); k++)
    {
        WaveFunctions.push_back(Shoot(EigenEnergies[k], ZStruct, Potential));
    }

    return(WaveFunctions);
}

