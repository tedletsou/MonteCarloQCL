#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "KGridSet.h"
#include "PoissonSolver.h"
#include "Constants.h"



KGridStruct CreateKSpaceGrid(int KSize, PoissonResult PResult, DeckDataStuct DeckData)
{
	//Enforce KSize to be an odd value, thus the middle of KSize will be zero, and KMagVec will run from 0 to Kmax
	if (KSize % 2 == 0)
	{
		KSize++;
	}

	// Find Maximum Electron momentum Kmax, 
	// where Kmax is found from Emax the bounds of the Potential, Emax=max(Potential)-min(Potential)

	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());

	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Calculate Effective mass at Emax, as an upper bound
	std::vector<double> msz = CalcEffectiveMass(Emax, PResult.NewZStruct);

	//Average over z, for the effective mass ms
	double ms = std::accumulate(msz.begin(), msz.end(), 0.0) / msz.size();

	//KMax is found using parabolic band
	double Kmax = sqrt(Emax * 2 * ms * ec) / hbar;

	//Initialize The KGrid Matricies
	std::vector<std::vector<double>> KGridKx(KSize, std::vector<double>(KSize, 0));
	std::vector<std::vector<double>> KGridKy(KSize, std::vector<double>(KSize, 0));
	std::vector<std::vector<double>> KGridKMag(KSize, std::vector<double>(KSize, 0));

	for (int n = 0; n < KSize; n++)
	{
		for (int m = 0; m < KSize; m++)
		{
			KGridKx[n][m] = double(m) / double(KSize-1) * 2 * Kmax - Kmax;
			KGridKy[n][m] = double(n) / double(KSize-1) * 2 * Kmax - Kmax;
			KGridKMag[n][m] = sqrt(KGridKx[n][m]* KGridKx[n][m]+ KGridKy[n][m] * KGridKy[n][m]);
		};
	};

	std::vector<double> KMagVec;

	// Find a vectorized magnitude of K, for ease of integration in other functions, i.e. |K| = |K|min : |K|max
	// This implementation avoids redundant values in the KmagVec
	for (int n = (KSize - 1) / 2; n < KSize; n++)
	{
		for (int m = n; m < KSize; m++)
		{
			KMagVec.push_back(sqrt(KGridKx[n][n] * KGridKx[n][n] + KGridKy[m][m] * KGridKy[m][m]));
		}
	}

	//Sort the Values of KMagVec from low to high
	std::sort(KMagVec.begin(), KMagVec.end());

	//For Debug
	/*
	std::cout << "KMagVec: " << std::endl;
	for (int n = 0; n < KMagVec.size(); n++)
	{
		std::cout << n << "  :  " << KMagVec[n] << std::endl;
	}
	std::cout << std::endl;
	std::cout << Kmax <<std::endl;

	std::cout << "KGridKMag: " << std::endl;
	for (int n = 0; n < KGridKMag.size(); n++)
	{
		for (int m = 0; m < KGridKMag.size(); m++)
		{
		std::cout << n << "," << m << "  :  " << KGridKMag[n][m] << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << Kmax << std::endl;
	*/
	KGridStruct KGrid;

	KGrid.KGridKx = KGridKx;
	KGrid.KGridKy = KGridKy;
	KGrid.KGridKMag = KGridKMag;
	KGrid.KMagVec = KMagVec;

		return KGrid;
}
