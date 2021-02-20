#ifndef CHARGEDENSITYCALC_H
#define CHARGEDENSITYCALC_H

#include <iostream>
#include <fstream>
#include <vector>
#include "GenerateZSpace.h"

//Struct containg vectors for the charge density and Fermi Level along Z
struct ChargeDistSturct 
{
	std::vector<double> ChargeDistZ;
	std::vector<double> FermiLevelZ;

};


// yet to be determined
ChargeDistSturct  CalcInitDopantDensity(ZMaterialParmsStruct ZStruct, double TL);


#endif