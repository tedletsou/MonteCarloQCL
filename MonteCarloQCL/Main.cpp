#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "BoundEigenDebug.h"
#include "ChargeDensityCalc.h"
#include "QCLMath.h"
#include "Shoot.h"
#include "BoundEigenenergies.h"
#include "EigenEnergyCalc.h"
#include "CalcWaveFunction.h"
#include "PoissonSolver.h"

int main() 
{
	//Parse the Input Data for the simulation from mcpp_input.dat
	DeckDataStuct DeckInput = Parse("mcpp_input.dat");

	//Create Vectors for the material parameters along the Z, Vectors stored in Struct ZMaterialParmsStruct
	ZMaterialParmsStruct ZMaterialStruct = CreateZParams(DeckInput);

	for (int AppEindexa =0; AppEindexa < DeckInput.field_vals.size(); AppEindexa++)
	{
		//Calculate the initial Potential of the QCL Structure with applied Bias and Conduction Band edge
		ZMaterialStruct = CalcPotential(ZMaterialStruct, DeckInput.field_vals[AppEindexa]);

		//!!!!!!!!! HARD CODED TEMP NEEDS to CHANGE
		double TL = 300;

		//Calculate initial Dopant Ion Distribution and Fermi Levels from Dopant Profile and Temperature 
		ChargeDistSturct IonizedDopantDensity = CalcInitDopantDensity(ZMaterialStruct, TL);
				
		// Find the Energy Bounds for each State in the Conduction Band
		QCLMat EnergyBounds = CalcEnergyBounds(ZMaterialStruct);

		
		//Allowed Error in Eigen Energies for Bound States
		double EnergyTolerance = 1e-8;

		// Find the Energy of each State in the Conduction Band based of Energy Bounds found above, using root finder in GNU Scientific Library (GSL), [Must inlcude in Project to function]
		std::vector<double> EigenEnergies = EigenEnergyCalc(EnergyBounds, ZMaterialStruct, EnergyTolerance);

		//Print out initial Eigen Energies
		std::cout << std::endl << "Initial Bound States  " << std::endl;
		for (int n = 0; n < EigenEnergies.size(); n++)
		{
			std::cout << EigenEnergies[n] << std::endl;
		}

		//Find Initial Wavefunctions based on Calculated Eigen Energies
		std::vector<WFStruct> WaveFunctions = CalculateWaveFunctions(EigenEnergies, ZMaterialStruct);

		//Calculate initial Fermilevel, Carrier Density in each subband and total charge density along Z (rho), returns only rho 
		std::vector<double> rho = CalcInitCarrierDensity(IonizedDopantDensity, ZMaterialStruct, WaveFunctions, EigenEnergies, TL);

		//Optional Code to write ZGrid and layer stucture to a File

		FILE* fpWF = fopen("WaveFunctions.txt", "w+");

		for (int n = 0; n < EigenEnergies.size(); n++)
		{
			fprintf(fpWF, "%.20g \t", EigenEnergies[n]);
		}

		fprintf(fpWF, "\n");

		for (int k = 0; k < WaveFunctions[0].Wavefunction.size(); k++)
		{
			for (int n = 0; n < WaveFunctions.size(); n++)
			{
				fprintf(fpWF, "%.20g \t", WaveFunctions[n].Wavefunction[k]);
			}
			fprintf(fpWF, "\n");
		}

		fclose(fpWF);


		//Optional Code to write ZGrid and layer stucture to a File

		FILE* fpCBE = fopen("Potential.txt", "w+");

		for (int k = 0; k < ZMaterialStruct.Potential.size(); k++)
		{
			fprintf(fpCBE, "%f \t", ZMaterialStruct.Potential[k]);
			fprintf(fpCBE, "%f \n", ZMaterialStruct.ZGridm[k]);
		}

		fclose(fpCBE);
			
	};



return 0;
}