#pragma once
#include <vector>
#include "GenerateZSpace.h"
#include "ParseInput.h"
#include "CalcWaveFunction.h"


std::vector<double> CalcEffectiveMass(double EnergyWF, ZMaterialParmsStruct ZStruct);

std::vector<double> CalcWeightedEffectiveMass(ZMaterialParmsStruct ZStruct, std::vector<WFStruct> Wavefunctions, std::vector<double> EigenEnergy);