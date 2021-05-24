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

double CalcAverageEffectiveMass(ZMaterialParmsStruct ZStruct, WFStruct Wavefunction, double Ei)
{
	//Calculates an effectice mass based on Non-Parabolicity, Follows "Modeling Techniques for Quantum Cascade Lasers" approach to non-parabolicty
	//Effective mass ms
	double msE=0;
	
	//Integrate over Z: this integral: ms*(1 + 2*alpha+beta)(Ei-V(z)) * Wavefunction(z) * Wavefunction(z)*DZ
	//Delta Z
	double DZ;

	//Independent variable of integration Z
	std::vector<double> Z = ZStruct.ZGridm;

	//Indexing over z
	for (int m = 0; m < ZStruct.ZGrid.size(); m++)
	{
		// Delta Z used for integration
		DZ = Z[m + 1] - Z[m];

		msE += ZStruct.ZMass[m] * (1 + (2 * ZStruct.ZAlphaNp[m] + ZStruct.ZBetaNp[m]) * (Ei - ZStruct.Potential[m])) * Wavefunction.Wavefunction[m] * Wavefunction.Wavefunction[m] * DZ;
	}

	return msE;
}


std::vector<double> CalcAllEffectiveMass(ZMaterialParmsStruct ZStruct, std::vector<WFStruct> Wavefunctions, std::vector<double> EigenEnergy)
{
	
	//Calculates the Effective Mass for the bottom of each subband, returns a vector of masses

	//Effective Mass vector m* initialized to 0
	std::vector<double> ms(Wavefunctions.size(), 0);

	//Integrate over Z: this integral: ms*(1 + 2*alpha+beta)(Ei-V(z)) * Wavefunction(z) * Wavefunction(z)*DZ
	//Delta Z
	double DZ;

	//Independent variable of integration Z
	std::vector<double> Z = ZStruct.ZGridm;

	for (int n = 0; n < Wavefunctions.size(); n++)
	{
		//Indexing over z
		for (int m = 0; m < ZStruct.ZGrid.size(); m++)
		{
			// Delta Z used for integration
			DZ = Z[m + 1] - Z[m];

			ms[n] += ZStruct.ZMass[m] * (1 + (2 * ZStruct.ZAlphaNp[m] + ZStruct.ZBetaNp[m]) * (EigenEnergy[n] - ZStruct.Potential[m])) * Wavefunctions[n].Wavefunction[m] * Wavefunctions[n].Wavefunction[m] * DZ;
		}
	}

	std::cout << std::endl;
	std::cout << "Effective Mass Printout"<< std::endl;
	for (int n = 0; n < EigenEnergy.size(); n++)
	{
		std::cout << "Band: " << n << "  ms: " << ms[n] << std::endl;
	}

	return(ms);
	
	//Original Function
	/*
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

	*/
}
