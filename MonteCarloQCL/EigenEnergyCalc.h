#pragma once

#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include "QCLMath.h"
#include "GenerateZSpace.h"


std::vector<double> EigenEnergyCalc(QCLMat EigenEnergyBounds, ZMaterialParmsStruct ZStruct, std::vector<double> Potential, double EnergyTolerance);

struct SolverParams 
{
	std::vector<double> Potential;
	ZMaterialParmsStruct ZStruct;
};
