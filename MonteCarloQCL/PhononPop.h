#pragma once

// Struct Containing the Parameters for Phonons, including occupancy number, Epsilon at Infinity, Epsilon static, LO phonon Energy
struct LOPhonStruct
{
	double Epinf;
	double Eps;
	double ELO;
	double NLO;
};

//Function that returns the parameters for LO Phonon Population in the LOPhonon Struct
LOPhonStruct LOPhonGaAsOccupancy(double TL);