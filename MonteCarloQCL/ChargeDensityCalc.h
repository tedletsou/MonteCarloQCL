#ifndef CHARGEDENSITYCALC_H
#define CHARGEDENSITYCALC_H

#include <iostream>
#include <fstream>
#include <vector>
#include "GenerateZSpace.h"
#include "CalcWaveFunction.h"

//Struct containg vectors for the charge density and Fermi Level along Z
struct ChargeDistSturct 
{
	std::vector<double> ChargeDistZ;
	std::vector<double> FermiLevelZ;

};

struct FermiSolverParams
{
	//EigenEnergies of the subbands
	std::vector<double> EigenEnergies;
	//Temperature of the electron distribution
	double T;
	//Effective mass for each subband
	std::vector<double> ms;
	//Total Carrier density given by charge neutrality with the Dopants Nd
	double Nstot;
};

// Function used by fsolver to calculate the fermi level, the function returs the sheet density for the subband Energies in params
double FermiLevel2DSheetDensity(double u, void* params);

// Calculated using Boltzman Distribution
ChargeDistSturct  CalcInitDopantDensity(ZMaterialParmsStruct ZStruct, double TL);

// Calculation using Thermal distibution of 2D sheet density with subband energies
std::vector<double> CalcInitCarrierDensity(ChargeDistSturct IonizedDopantDensity, ZMaterialParmsStruct ZStruct, std::vector<WFStruct> WaveFunctions, std::vector<double> EigenEnergies, double TL);

// Function used by fsolver to calculate the fermi level
double CalcFermiLevel(ChargeDistSturct IonizedDopantDensity, ZMaterialParmsStruct ZStruct, std::vector<WFStruct> WaveFunctions, std::vector<double> EigenEnergies, double TL);


#endif