#pragma once

#include "Shoot.h"
#include "QCLMath.h"

std::vector<WFStruct> CalculateWaveFunctions(std::vector<double> EigenEnergies, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms);
