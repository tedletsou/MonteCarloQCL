#include <vector>
#include "EffectiveMass.h"
#include "GenerateZSpace.h"
#include "ParseInput.h"
#include "Constants.h"
#include <math.h>
#include <iostream>

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