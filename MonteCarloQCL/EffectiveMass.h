#pragma once
#include <vector>
#include "GenerateZSpace.h"
#include "ParseInput.h"
#include "CalcWaveFunction.h"


std::vector<double> CalcEffectiveMass(double EnergyWF, ZMaterialParmsStruct ZStruct);

std::vector<double> CalcAllEffectiveMass(ZMaterialParmsStruct ZStruct, std::vector<WFStruct> Wavefunctions, std::vector<double> EigenEnergy);

//Can only be used after calculating potential and Wavefunction
double CalcAverageEffectiveMass(ZMaterialParmsStruct ZStruct,WFStruct Wavefunction, double Ei);