#include <vector>
#include "EffectiveMass.h"
#include "GenerateZSpace.h"
#include "ParseInput.h"
#include "Constants.h"
#include <math.h>
#include <iostream>
#include "Shoot.h"

std::vector<double> CalcEffectiveMass(double EnergyWF,ZMaterialParmsStruct ZStruct) 
{
	//Effective Mass vector m*
	std::vector<double> Ms{};

	//Iterate through Z grid creating effective mass for each point in Z
	for (int k = 0; k < ZStruct.ZMass.size(); k++)
	{
		//Ms.push_back(me * 1 / ( (2 / 3) * ZStruct.ZKane[k] / (EnergyWF - ZStruct.ZLHole[k]) + (1 / 3) * ZStruct.ZKane[k] / (EnergyWF - ZStruct.ZSploff[k])));
		Ms.push_back(me / ((0.6666666)*ZStruct.ZKane[k] / (EnergyWF-ZStruct.ZLHole[k]) + (.33333333) * ZStruct.ZKane[k] / (EnergyWF - ZStruct.ZSploff[k])));
	}

	return(Ms);
}

std::vector<double> CalcWeightedEffectiveMass(ZMaterialParmsStruct ZStruct, std::vector<WFStruct> Wavefunctions, std::vector<double> EigenEnergy)
{
	//Effective Mass vector m* initialized to 0
	std::vector<double> ms(EigenEnergy.size(),0);

	//Implement the following integral: ms * {1 + [2*alpha + Beta] [Ei-V] Abs(Wavefunction(x)}

	//Iterate through each EigenEngergy
	for (int n = 0; n < EigenEnergy.size(); n++)
	{
		//Iterate through Z grid
		for (int k = 0; k < ZStruct.ZMass.size(); k++)
		{
			// 3 Band k dot p effective mass model using Kane Energy
			ms[n] += me / ((0.6666666) * ZStruct.ZKane[k] / (EigenEnergy[n] - ZStruct.ZLHole[k]) + (.33333333) * ZStruct.ZKane[k] / (EigenEnergy[n] - ZStruct.ZSploff[k]));

			//Modify ms based on the Potential and Wavefunctions
			//ms[n] = Wavefunctions[n].Wavefunction[k] * Wavefunctions[n].Wavefunction[k] * (EigenEnergy[n] - ZStruct.Potential[k]);
		}
	}

	return(ms);
}