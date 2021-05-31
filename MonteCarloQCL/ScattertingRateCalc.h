#pragma once

#include <vector>
#include "FormFactorCalc.h"
#include "PoissonSolver.h"
#include "PhononPop.h"
#include "KGridSet.h"

typedef std::vector<std::vector<std::vector<double>>> ScatteringRateMatrix;

typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ScatteringRateEEMatrix;

ScatteringRateMatrix LOPhononEmitScatRateCalc(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL);

ScatteringRateMatrix LOPhononAbsScatRateCalc(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL);

ScatteringRateEEMatrix EEScatRateCalc(FormFactorEEStruct EEFF, PoissonResult PResult, KGridStruct KGrid, double TL, int Numq);

