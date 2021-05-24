#pragma once

#include <vector>
#include "PoissonSolver.h"
#include "PhononPop.h"

typedef std::vector<std::vector<std::vector<double>>> FormFactorMatrix;

typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> FormFactorEEMatrix;

struct FormFactorStruct
{
	FormFactorMatrix FormFactor;
	std::vector<double> qvals;
};

struct FormFactorEEStruct
{
	FormFactorEEMatrix FormFactor;
	std::vector<double> qvals;
};

double FormFactorInterpolate(FormFactorStruct FF, int i, int f, double q);

double FormFactorEEInterpolate(FormFactorEEStruct FF, int i, int f, int g, int H, double qxy);

FormFactorEEStruct FormFactorUpSample(FormFactorEEStruct FF, int UpFactor);

FormFactorStruct FormFactorLOPhononCalc(PoissonResult PResult, LOPhonStruct PhononParam, int Numq);

FormFactorEEStruct FormFactorEECalc(PoissonResult PResult, int Numq);