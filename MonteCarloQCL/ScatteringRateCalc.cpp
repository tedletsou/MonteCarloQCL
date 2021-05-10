#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "EffectiveMass.h"
#include "FormFactorCalc.h"
#include "PoissonSolver.h"
#include "Constants.h"
#include "PhononPop.h"
#include "KGridSet.h"
#include "ScattertingRateCalc.h"

ScatteringRateMatrix LOPhononEmitScatRate(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL)
{
	//Calculate LO Phonon Parameters for GaAs
	LOPhonStruct LOGaAs = LOPhonGaAsOccupancy(TL);

	//Number of Subbands, used for convenience
	int NumSb = PResult.NewWaveFunctions.size();

	//Number of initial k-values to iterate over
	int Numk = KGrid.KMagVec.size();

	//Number of Theta values to integrate over
	int NumTh = LOPhononFF.qvals.size();

	//Delta Theta 
	double DThetha = 2 * Pi / (NumTh-1);
	
	//independent variable used for integration from0 to 2*Pi
	double Theta;

	//Momentum of LO Phonon involved in scattering
	double q;

	//Final Kf^2 of scattered electron
	double Kf;

	//Compute Scattering Rate Integral for initial subband i, final subband f, and magnitude of initial electron momentum ki
	
	//Initialize Scattering Rate Matrix zero
	ScatteringRateMatrix Result(NumSb, std::vector<std::vector<double>>(NumSb, std::vector<double>(Numk, 0)));

	//Compute Effective Mass
	std::vector<double> msz = CalcEffectiveMass(0, PResult.NewZStruct);
	double ms = *std::min_element(msz.begin(), msz.end());

	//Initialize Interpolated FormFactor Value
	double FFInterp;

	//Indexing over initial subband i
	for (int i = 0; i < NumSb; i++)
	{
		//Indexing over final subband f
		for (int f = 0; f < NumSb; f++)
		{
			//Indexing over magnitude of ki	
			for (int n = 0; n < Numk; n++)
			{
				//Initial energy of electron for non-parabolicity 
				double Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms*ec);
				
				//Calculate Effective mass to account for non-parabolicity for given Ei
				msz = CalcEffectiveMass(Ei, PResult.NewZStruct);
				ms = *std::min_element(msz.begin(), msz.end());

				//Find initial energy of electron accounting for non-parabolicity, i.e. using the proper ms 
				Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms * ec);

				if (Ei - PResult.EigenEnergies[f] > LOGaAs.ELO)
				{
					for (int th = 0; th < NumTh; th++)
					{

						//Theta the independent integration variable
						Theta = th * DThetha;

						//Find Kf the magnitude of the final K vector using: sqrt( ki^2 + (2ms)/h^2 * (Ei - Ef -ELO) )
						Kf = sqrt( KGrid.KMagVec[n] * KGrid.KMagVec[n] + 2 * ms *ec / (hbar * hbar) * (PResult.EigenEnergies[i] - PResult.EigenEnergies[f] - LOGaAs.ELO));

						//Find the magnitude of q the phonon momentum using: sqrt( ki^2 + kf^2 - 2*ki*kf*cos(Theta) )
						q = sqrt(KGrid.KMagVec[n] * KGrid.KMagVec[n] + Kf * Kf - 2 * KGrid.KMagVec[n] * Kf * cos(Theta));

						//Using the q calculated above, interpolate the FormFactor for LO Phonon Scattering
						FFInterp = FormFactorInterpolate(LOPhononFF, i, f, q);

						//Accumulate the Result of: ms*e^2*ELO/(8*pi*hbar^3)(1/Epinf - 1/Eps)*(NLO +1) integral( A(q)/q dTHeta from 0 to 2*pi)
						Result[i][f][n] += ms * ec * ec * LOGaAs.ELO * ec / (8 * Pi * hbar * hbar * hbar * Ep0) *(1 / LOGaAs.Epinf - 1 / LOGaAs.Eps)* (LOGaAs.NLO + 1) *FFInterp / q * DThetha;
						
					}
					//Printing Vals for Debug
					//std::cout << "i: " << i << "  f: " << f << "  MaxK: " << KGrid.KMagVec[Numk-1] << "  Ki: " << KGrid.KMagVec[n] << "  q: " << q << "  Ei: " << PResult.EigenEnergies[i] << "  Ef: " <<PResult.EigenEnergies[f] << "  ms: " << ms << "  FFInterp: " << FFInterp << "  Result: " << Result[i][f][n] << std::endl;
				}
				else
				{
					Result[i][f][n] = 0;
				}
				
			}
		}
	}

	return Result;
}