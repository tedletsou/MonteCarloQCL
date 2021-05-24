#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "PoissonSolver.h"
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
		ZStruct.Potential[k] = ZStruct.CBand[k] + NewPotential[k];
	}
	
	ZMaterialParmsStruct ZStructPot = ZStruct;

	/*
	//File Output of the Efield from Poisson
	FILE* fpPE = fopen("EFieldPoiss.txt", "w+");

	for (int k = 0; k < EField.size(); k++)
	{
		fprintf(fpPB, "%f \t", EField[k]);
		fprintf(fpPB, "%f \n", ZStruct.ZGridm[k] * 1e10);
	}

	fclose(fpPB);
	*/

	return(ZStructPot);
};


PoissonResult PoissonSolver(ZMaterialParmsStruct OldZStruct, ChargeDistSturct IonizedDopantDensity, std::vector<double> InitRho, std::vector<WFStruct> InitWaveFunctions, double TL, double ErrorTol, double AppField)
{
	int counter = 0;
	double Error=0;
	
	//FermiLevel
	double u;
	
	struct PoissonResult Result;

	//Use Initial calculation of rho in Zsturt
	OldZStruct.rho = InitRho;

	//Initialize New Zstruct and Wavefunctions
	ZMaterialParmsStruct NewZStruct= OldZStruct;
	std::vector<WFStruct> NewWaveFunctions;

	do
	{
		//Increment Counter
		counter++;
		
		//Calculate the Potential of the QCL Structure with applied Bias and Conduction Band edge with the initial calculation of rho
		NewZStruct = CalcPotential(NewZStruct, AppField);

		// Find the New Energy Bounds for the new Potential
		QCLMat NewEnergyBounds = CalcEnergyBounds(NewZStruct);

		//Allowed Error in Eigen Energies for Bound States
		double EnergyTolerance = 1e-8;

		// Find the Energy of each State in the Conduction Band based of Energy Bounds found above, using root finder in GNU Scientific Library (GSL), [Must inlcude in Project to function]
		std::vector<double> NewEigenEnergies = EigenEnergyCalc(NewEnergyBounds, NewZStruct, EnergyTolerance);

		//Print out Eigen Energies
		std::cout << std::endl << "Initial Bound States  " << std::endl;
		for (int n = 0; n < NewEigenEnergies.size(); n++)
		{
			std::cout << NewEigenEnergies[n] << std::endl;
		}

		//Find Initial Wavefunctions based on Calculated Eigen Energies
		std::vector<WFStruct> NewWaveFunctions = CalculateWaveFunctions(NewEigenEnergies, NewZStruct);

		//Update the Old Zstruct with the New ZStruct
		OldZStruct.rho = NewZStruct.rho;

		//Calculate an update for charge density based on New Wavefunctions
		NewZStruct.rho = CalcInitCarrierDensity(IonizedDopantDensity, NewZStruct, NewWaveFunctions, NewEigenEnergies, TL);
		u = CalcFermiLevel(IonizedDopantDensity, NewZStruct, NewWaveFunctions, NewEigenEnergies, TL);

		//Calculate the RMS difference of the Old and New values of rho in OldZStruct and NewZStruct along Z
		for (int k = 0; k < NewZStruct.rho.size(); k++)
		{
			Error += (NewZStruct.rho[k] - OldZStruct.rho[k]) * (NewZStruct.rho[k] - OldZStruct.rho[k]);
		}
		Error = sqrt(Error / (NewZStruct.rho.size()));

		//For Debug write Counter
		std::cout << "Counter: " << counter << "  Error: " << Error << std::endl;

		Result.NewWaveFunctions = NewWaveFunctions;
		Result.NewZStruct = NewZStruct;
		Result.EigenEnergies = NewEigenEnergies;
		Result.u = u;
	}
	while (Error > ErrorTol);
	
	return(Result);
	
};