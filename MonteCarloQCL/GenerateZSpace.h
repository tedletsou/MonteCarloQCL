#ifndef GENERATEZSPACE_H
#define GENERATEZSPACE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

//Thickness of Each layer along Z
vector<double> ZThickness;
//Intereface Posistions along the Z axis
vector<double> ZSpace;
//Grid along the growth direction used to define potentials, Effective Masses and Wavefunctions
vector <double> ZGrid;
//Conduction Band Edge CBE
vector<double> CBE;
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

int MeshDen;

#endif