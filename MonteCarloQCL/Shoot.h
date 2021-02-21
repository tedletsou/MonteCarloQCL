#pragma once

// Header file for Shoot, used to calculate Psi for an arbitrary potential, 

// Returns vector containg Wavefunction (Wf), double last value of Wf, int number of zeros in Wf

// Pass in Energy Value for Wf (WfEndegy), Potential, Material Paramets along Z axis contained in Struct ZMaterialParmsStruct

#include <vector>
#include "GenerateZSpace.h"


struct WFStruct 
{
	std::vector<double> Wavefunction;
	double LastValue;
	int NumZeros;
};

WFStruct Shoot(double WfEnergy, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms);