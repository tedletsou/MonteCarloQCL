#pragma once

#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include "QCLMath.h"
#include "GenerateZSpace.h"


std::vector<double> EigenEnergyCalc(QCLMat EigenEnergyBounds, ZMaterialParmsStruct ZStruct, double EnergyTolerance);

struct SolverParams 
{
	ZMaterialParmsStruct ZStruct;
};
