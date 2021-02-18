#include <iostream>
#include <fstream>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include <vector>

std::vector<double> CalcChargeDensity(ZMaterialParmsStruct ZStruct, double TL)
{
    //Effective Electron Density following Baird
    std::vector<double> G;

    G = MultiplyVectorByScalar(ZStruct.ZDoping, 0.5 * pow((2 * Pi * pow(hbar, 2) / (kb * TL)), 3 / 2));

    return std::vector<double>();
}
