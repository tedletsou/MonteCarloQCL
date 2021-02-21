#include <iostream>
#include <fstream>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include "Shoot.h"
#include <vector>
#include <cmath>
#include <algorithm>


double CalcEnergyBounds(ZMaterialParmsStruct ZStruct)
{
	// Find min and max conduction band energy 
	auto MinMax = std::minmax_element(ZStruct.CBand.begin(), ZStruct.CBand.end());
	
	// Calculate number of zeros at the min and max conduction band energy
	WFStruct BottomBand = Shoot(*MinMax.first, ZStruct.ZGridm, ZStruct.CBand, ZStruct.ZMass);
	WFStruct TopBand = Shoot(*MinMax.second, ZStruct.ZGridm, ZStruct.CBand, ZStruct.ZMass);

	// The number of wavefunctions is equal to the number of zeros at the top of the band - numbe of zeros at bottom
	int NumWavefunction = TopBand.NumZeros - BottomBand.NumZeros;

	// Initializing matrix that keeps track of zeros
	QCLMat ZeroAtBounds{ std::vector<double> (NumWavefunction, 0), std::vector<double>(NumWavefunction, NumWavefunction) };

	// Initializing matrix that keeps track of bounded energies
	QCLMat EnergyAtBounds{ std::vector<double>(NumWavefunction, 0), std::vector<double>(NumWavefunction, *MinMax.second) };

	// Initial midpoint energy for bisection search of eigenenergy bounds
	double MidpointEnergy = (*MinMax.second - *MinMax.first) / 2;

	for (int k = 0; k < NumWavefunction - 1; k++) 
	{
		while (ZeroAtBounds[1][k] - ZeroAtBounds[0][k] > 1)
		{
			WFStruct MidPoint = Shoot(MidpointEnergy, ZStruct.ZGridm, ZStruct.CBand, ZStruct.ZMass);
			std::cout << MidPoint.NumZeros;

		}
	}
	return 0;

}