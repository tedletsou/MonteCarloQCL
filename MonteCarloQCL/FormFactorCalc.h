#pragma once

#include <vector>
#include "PoissonSolver.h"
#include "PhononPop.h"

typedef std::vector<std::vector<std::vector<double>>> FormFactorMatrix;

struct FormFactorStruct
{
	FormFactorMatrix FormFactor;
	std::vector<double> qvals;
};

double FormFactorInterpolate(FormFactorStruct FF, int i, int f, double q);

FormFactorStruct FormFactorLOPhononCalc(PoissonResult PResult, LOPhonStruct PhononParam, int Numq);