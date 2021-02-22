#pragma once

#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include "QCLMath.h"


std::vector<double> EigenEnergyCalc(QCLMat EigenEnergyBounds, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms, double EnergyTolerance);

struct SolverParams 
{
	std::vector<double> ZGridm;
	std::vector<double> Potential;
	std::vector<double> ms;
};
