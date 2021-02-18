#ifndef GENERATEZSPACE_H
#define GENERATEZSPACE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "GenerateZSpace.h"
#include "ParseInput.h"

using namespace std;


struct ZMaterialParmsStruct {
	vector <double> ZGrid;
	//Conduction Band Edge CBE
	vector<double> CBand;
	//Effective Masses along Z
	vector<double> ZMass;
	//Valence Band Edge VBand along Z
	vector<double> ZVBand;
	//Light Hole band Edge VBand along Z
	vector<double> ZLHole;
	//Split off band Edge VBand along Z
	vector<double> ZSploff;
	// Delta Energy Difference between Light Hole and Split Off Bands along Z
	vector<double> ZDelta;
	// Kane Energies along Z
	vector<double> ZKane;
	// Band Gap along Z
	vector<double> ZBandGap;
	// Doping along Z
	vector<double> ZDoping;
};

ZMaterialParmsStruct CreateZParams(DeckDataStuct DeckData);


#endif