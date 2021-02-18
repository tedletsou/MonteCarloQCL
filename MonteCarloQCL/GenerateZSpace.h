#ifndef GENERATEZSPACE_H
#define GENERATEZSPACE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "GenerateZSpace.h"
#include "ParseInput.h"

// Define output struct for CreateZParams
struct ZMaterialParmsStruct {
	std::vector <double> ZGrid;
	//Conduction Band Edge CBE
	std::vector<double> CBand;
	//Effective Masses along Z
	std::vector<double> ZMass;
	//Valence Band Edge VBand along Z
	std::vector<double> ZVBand;
	//Light Hole band Edge VBand along Z
	std::vector<double> ZLHole;
	//Split off band Edge VBand along Z
	std::vector<double> ZSploff;
	// Delta Energy Difference between Light Hole and Split Off Bands along Z
	std::vector<double> ZDelta;
	// Kane Energies along Z
	std::vector<double> ZKane;
	// Band Gap along Z
	std::vector<double> ZBandGap;
	// Doping along Z
	std::vector<double> ZDoping;
};

// Define create zparams function
ZMaterialParmsStruct CreateZParams(DeckDataStuct DeckData);


#endif