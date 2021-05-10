#include <cmath>
#include "PhononPop.h"
#include "Constants.h"

LOPhonStruct LOPhonGaAsOccupancy(double TL)
{
	// Struct Containing the Parameters for GaAs Phonon, including occupancy number
	LOPhonStruct GaAsLOPhon;

	// LO Phonon Energy at Gamma Point for GaAs
	GaAsLOPhon.ELO = .036;

	//Static and High Frequency Permitivity from www.offe.ru/SVA/NSM/Semicond/GaAs/basic.html
	GaAsLOPhon.Epinf = 10.89;
	GaAsLOPhon.Eps = 12.9;

	//Occupancy of the Phonon LO phonon mode, given by Bose-Einstien Distribution
	GaAsLOPhon.NLO = 1 / exp(GaAsLOPhon.ELO / (kb * TL));

	return GaAsLOPhon;
}