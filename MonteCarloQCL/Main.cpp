#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "ChargeDensityCalc.h"
#include "QCLMath.h"
#include <fstream>

int main() 
{
	//Parse the Input Data for the simulation from mcpp_input.dat
	DeckDataStuct DeckInput = Parse("mcpp_input.dat");

	//Create Vectors for the material parameters along the Z, Vectors stored in Struct ZMaterialParmsStruct
	ZMaterialParmsStruct ZMaterialStruct = CreateZParams(DeckInput);

	double TL = 10;

	CalcChargeDensity(ZMaterialStruct, TL);

return 0;
}