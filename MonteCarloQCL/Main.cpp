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

int main() 
{
	//Parse the Input Data for the simulation from mcpp_input.dat
	DeckDataStuct DeckInput = Parse("mcpp_input.dat");

	//Create Vectors for the material parameters along the Z, Vectors stored in Struct ZMaterialParmsStruct
	ZMaterialParmsStruct ZMaterialStruct = CreateZParams(DeckInput);

	double TL = 10;

	//ZMaterialParmsStruct ZMaterialStruct2 = ZMaterialStruct;

	
	//Calculate initial Dopant Ion Distribution and Fermi Levels from Dopant Profile and Temperature 
	ChargeDistSturct InitDopantDensity = CalcInitDopantDensity(ZMaterialStruct, TL);

	
	//Calculate Band Structure
	double E = 0.5;
	
	WFStruct Wf = Shoot(E, ZMaterialStruct.ZGridm, ZMaterialStruct.CBand, ZMaterialStruct.ZMass);

	

	for (int k = 0; k < Wf.Wavefunction.size(); k++)
	{
		std::cout << Wf.Wavefunction[k] << std::endl;
	}

	std::cout << Wf.NumZeros << std::endl;

	

return 0;
}