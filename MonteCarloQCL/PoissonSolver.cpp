#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include "Shoot.h"
#include "BoundEigenenergies.h"
#include "EigenEnergyCalc.h"


ZMaterialParmsStruct CalcPotential(ZMaterialParmsStruct ZStruct, double AppliedField)
{
	//Function Updates Zstruct with a New potential based on Poisson Effect, Applied Voltage and Conduction Band Edge

	std::vector<double> EField;
	EField.resize(ZStruct.ZGrid.size());

	std::vector<double> NewPotential;
	NewPotential.resize(ZStruct.ZGrid.size());

	double DZ;
	
	//Integrate rho/Permitivity from 0 to z for each value along the z grid to get Electrict Field along Z grid
	for (int k = 0; k < ZStruct.ZGrid.size(); k++)
	{
		for (int m = 0; m <= k - 1; m++)
		{
			//Delta Z used for integration
			DZ = ZStruct.ZGridm[m + 1] - ZStruct.ZGridm[m];

			EField[k] += ZStruct.rho[m]/(ZStruct.ZPermitivity[m]* Ep0) * DZ;
		}
		//
		EField[k] += AppliedField;
	}

	//Integrate Electric Field from 0 to z for each value along the z grid to get Potenital along Z grid
	for (int k = 0; k < ZStruct.ZGrid.size(); k++)
	{
		for (int m = 0; m <= k - 1; m++)
		{
			//Delta Z used for integration
			DZ = ZStruct.ZGridm[m + 1] - ZStruct.ZGridm[m];

			//In untis of Volts
			NewPotential[k] += EField[m]*DZ;
		}
	}

	//Add the Potential from rho and the applied field to the Conduction Band edge
		//Integrate Electric Field from 0 to z for each value along the z grid
	for (int k = 0; k < ZStruct.ZGrid.size(); k++)
	{
		ZStruct.Potential[k] += NewPotential[k];
	}
	
	ZMaterialParmsStruct ZStructPot = ZStruct;

	return(ZStructPot);
};


ZMaterialParmsStruct PoissonSolver(ZMaterialParmsStruct ZStruct, ChargeDistSturct IonizedDopantDensity, double TL)
{
	/*
	// Find the Energy Bounds for each State in the Conduction Band
	QCLMat EnergyBounds = CalcEnergyBounds(ZStruct);

	//Allowed Error in Eigen Energies for Bound States
	double EnergyTolerance = 1e-8;

	// Find the Energy of each State in the Conduction Band based of Energy Bounds found above, using root finder in GNU Scientific Library (GSL), [Must inlcude in Project to function]
	std::vector<double> EigenEnergies = EigenEnergyCalc(EnergyBounds, ZStruct, EnergyTolerance);

	//Print out initial Eigen Energies
	std::cout << std::endl << "Initial Bound States  " << std::endl;
	for (int n = 0; n < EigenEnergies.size(); n++)
	{
		std::cout << EigenEnergies[n] << std::endl;
	}

	//Find Initial Wavefunctions based on Calculated Eigen Energies
	std::vector<WFStruct> WaveFunctions = CalculateWaveFunctions(EigenEnergies, ZMaterialStruct, ZMaterialStruct.CBand);

	//Calculate initial Fermilevel, Carrier Density in each subband and total charge density along Z (rho), returns only rho
	std::vector<double> rho = CalcInitCarrierDensity(IonizedDopantDensity, ZMaterialStruct, WaveFunctions, EigenEnergies, TL);

	*/
};