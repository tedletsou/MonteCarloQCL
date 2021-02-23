#include <iostream>
#include <fstream>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include "Shoot.h"
#include <vector>
#include <cmath>
#include <algorithm>


QCLMat CalcEnergyBounds(ZMaterialParmsStruct ZStruct)
{
	//Used to find the Energy Bounds for each state in the Conduction Band, uses number of Zeros in the Wavefunction to distinguish between states

	// Find min and max conduction band energy 
	auto MinMax = std::minmax_element(ZStruct.CBand.begin(), ZStruct.CBand.end());
	
	// Calculate number of zeros at the min and max conduction band energy
	WFStruct BottomBand = Shoot(*MinMax.first, ZStruct, ZStruct.CBand);
	WFStruct TopBand = Shoot(*MinMax.second, ZStruct, ZStruct.CBand);

	// The number of wavefunctions is equal to the number of zeros at the top of the band - numbe of zeros at bottom
	int NumWavefunction = TopBand.NumZeros - BottomBand.NumZeros;

	// Initializing matrix that keeps track of zeros
	QCLMat ZeroAtBounds{ std::vector<double> (NumWavefunction, 0), std::vector<double>(NumWavefunction, NumWavefunction) };

	// Initializing matrix that keeps track of bounded energies
	QCLMat EnergyAtBounds{ std::vector<double>(NumWavefunction, 0), std::vector<double>(NumWavefunction, *MinMax.second) };

	// Initial midpoint energy for bisection search of eigenenergy bounds
	double MidpointEnergy = (*MinMax.second + *MinMax.first) / 2;

	for (int k = 0; k < NumWavefunction; k++) 
	{
		while (ZeroAtBounds[1][k] - ZeroAtBounds[0][k] > 1)
		{
			// Generate Wave function at current Energy Bound
			WFStruct MidPointWF = Shoot(MidpointEnergy, ZStruct, ZStruct.CBand);

			//Interate through the EnergyBound and ZeroBounds Matricies and update with new bounds based on Energy from last shooting calculation
			for (int n = 0; n < NumWavefunction; n++)
			{
				//Check to see if the current state has less zeros than the number of Zeros of the Newly Calculated Wave Function & THe current State has MORE Energy than the Newly Calculated WF
				if ((n < MidPointWF.NumZeros) & (EnergyAtBounds[1][n] > MidpointEnergy))
				{
					EnergyAtBounds[1][n] = MidpointEnergy;
					ZeroAtBounds[1][n] = MidPointWF.NumZeros;
				}
				//Check to see if the current state has more zeros than the number of Zeros of the Newly Calculated Wave Function & THe current State has LESS Energy than the Newly Calculated WF
				if ((n >= MidPointWF.NumZeros) & (EnergyAtBounds[0][n] < MidpointEnergy))
				{
					EnergyAtBounds[0][n] = MidpointEnergy;
					ZeroAtBounds[0][n] = MidPointWF.NumZeros;
				}
			}

			//Calculate new Bisected Energy from previous bound
			MidpointEnergy = (EnergyAtBounds[1][k] + EnergyAtBounds[0][k]) / 2;
		
			std::cout << MidPointWF.NumZeros;

		}
	}

	//// For Test Printing of Bounds to States
	//for (int n = 0; n < 9; n++)
	//{
	//	std::cout  << std::endl << ZeroAtBounds[0][n] << std::endl;
	//}

	return EnergyAtBounds;

}