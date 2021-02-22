#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "ChargeDensityCalc.h"
#include "QCLMath.h"
#include "Shoot.h"
#include <fstream>
#include "BoundEigenenergies.h"
#include "EigenEnergyCalc.h"
#include "CalcWaveFunction.h"

int main() 
{
	//Parse the Input Data for the simulation from mcpp_input.dat
	DeckDataStuct DeckInput = Parse("mcpp_input.dat");

	//Create Vectors for the material parameters along the Z, Vectors stored in Struct ZMaterialParmsStruct
	ZMaterialParmsStruct ZMaterialStruct = CreateZParams(DeckInput);

	double TL = 10;

	//Calculate initial Dopant Ion Distribution and Fermi Levels from Dopant Profile and Temperature 
	ChargeDistSturct InitDopantDensity = CalcInitDopantDensity(ZMaterialStruct, TL);

	//Calculate Band Structure
	double E = 0.4377;
	
	// Find the Energy Bounds for each State in the Conduction Band
	QCLMat EnergyBounds = CalcEnergyBounds(ZMaterialStruct);

	//Allowed Error in Eigen Energies for Bound States
	double EnergyTolerance = 1e-8;


	std::cout << std::endl << "Energy Bounds" << std::endl;

	for (int n = 0; n < EnergyBounds[0].size(); n++)
	{
		std::cout << EnergyBounds[0][n] << EnergyBounds[1][n] << std::endl << std::endl;
	}


	// Find the Energy of each State in the Conduction Band based of Energy Bounds found above, using root finder in GNU Scientific Library (GSL), [Must inlcude in Project to function]
	std::vector<double> EigenEnergies = EigenEnergyCalc(EnergyBounds, ZMaterialStruct.ZGridm, ZMaterialStruct.CBand, ZMaterialStruct.ZMass, EnergyTolerance);

	std::cout << std::endl << "Energy Test Part" << std::endl;

	for (int n = 0; n < EigenEnergies.size(); n++)
	{
		std::cout << EigenEnergies[n] << std::endl;
	}

	std::vector<WFStruct> WaveFunctions = CalculateWaveFunctions(EigenEnergies, ZMaterialStruct.ZGridm, ZMaterialStruct.CBand, ZMaterialStruct.ZMass);


	//Optional Code to write ZGrid and layer stucture to a File
	/**/
	FILE* fpCBE = fopen("WaveFunctions.txt", "w+");

	for (int n = 0; n < EigenEnergies.size(); n++)
	{
		fprintf(fpCBE, "%.20g \t", EigenEnergies[n]);
	}

	fprintf(fpCBE, "\n");

	for (int k = 0; k < WaveFunctions[0].Wavefunction.size(); k++)
	{
		for (int n = 0; n < WaveFunctions.size(); n++)
		{
			fprintf(fpCBE, "%.20g \t", WaveFunctions[n].Wavefunction[k]);
		}
		fprintf(fpCBE, "\n");
	}

	fclose(fpCBE);

return 0;
}