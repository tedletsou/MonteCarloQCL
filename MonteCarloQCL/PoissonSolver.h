#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "GenerateZSpace.h"
#include "CalcWaveFunction.h"
#include "BoundEigenenergies.h"
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "EffectiveMass.h"
#include "EigenEnergyCalc.h"
#include "QCLMath.h"
#include "Shoot.h"


struct PoissonResult {
	std::vector<WFStruct> NewWaveFunctions;
	ZMaterialParmsStruct NewZStruct;
	std::vector<double> EigenEnergies;
};

ZMaterialParmsStruct CalcPotential(ZMaterialParmsStruct ZStruct, double AppliedField);

PoissonResult PoissonSolver(ZMaterialParmsStruct OldZStruct, ChargeDistSturct IonizedDopantDensity, std::vector<double> InitRho, std::vector<WFStruct> InitWaveFunctions, double TL, double ErrorTol, double AppField);