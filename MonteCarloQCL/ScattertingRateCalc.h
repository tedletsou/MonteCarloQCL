#pragma once

#include <vector>
#include "FormFactorCalc.h"
#include "PoissonSolver.h"
#include "PhononPop.h"
#include "KGridSet.h"

typedef std::vector<std::vector<std::vector<double>>> ScatteringRateMatrix;

ScatteringRateMatrix LOPhononEmitScatRate(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL);


