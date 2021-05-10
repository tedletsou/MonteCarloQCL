#pragma once
#include <vector>
#include "PoissonSolver.h"
#include "ParseInput.h"

struct KGridStruct
{
	std::vector<std::vector<double>> KGridKx;
	std::vector<std::vector<double>> KGridKy;
	std::vector<std::vector<double>> KGridKMag;
	std::vector<double> KMagVec;
};


KGridStruct CreateKSpaceGrid(int KSize, PoissonResult PResult, DeckDataStuct DeckData);