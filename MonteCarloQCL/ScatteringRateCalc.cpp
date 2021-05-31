#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "EffectiveMass.h"
#include "FormFactorCalc.h"
#include "PoissonSolver.h"
#include "Constants.h"
#include "PhononPop.h"
#include "KGridSet.h"
#include "QCLMath.h"
#include "ScattertingRateCalc.h"

ScatteringRateMatrix LOPhononEmitScatRateCalc(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL)
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

	//Find Maximum energy of electron before leaving the well, Emax
	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());

	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Compute Scattering Rate Integral for initial subband i, final subband f, and magnitude of initial electron momentum ki
	
	//Initialize Phonon Emission Scattering Rate Matrix zero
	ScatteringRateMatrix EmitResult(NumSb, std::vector<std::vector<double>>(NumSb, std::vector<double>(Numk, 0)));

	//Initialize Ei and ms, Initial Energy and Effective mass
	double Ei, ms;

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
				//Compute initial guess for Effective Mass
				ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions[i], PResult.EigenEnergies[i]);

				//Initial energy of electron for non-parabolicity 
				Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms*ec);
				
				//Calculate Effective mass to account for non-parabolicity for given Ei
				ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions[i], Ei);

				//Find initial energy of electron accounting for non-parabolicity, i.e. using the proper ms 
				Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms * ec);

				if ((Ei - PResult.EigenEnergies[f] >= LOGaAs.ELO) & (Ei <= Emax))
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
						EmitResult[i][f][n] += ms * ec * ec * LOGaAs.ELO * ec / (8 * Pi * hbar * hbar * hbar * Ep0) *(1 / LOGaAs.Epinf - 1 / LOGaAs.Eps)* (LOGaAs.NLO + 1) *FFInterp / q * DThetha;
						
					}
					//Printing Vals for Debug
					//std::cout << "NLO:  " << LOGaAs.NLO << std::endl;
					//std::cout << "i: " << i << "  f: " << f << "  MaxK: " << KGrid.KMagVec[Numk-1] << "  Ki: " << KGrid.KMagVec[n] << "  q: " << q << "  Ei: " << PResult.EigenEnergies[i] << "  Ef: " <<PResult.EigenEnergies[f] << "  ms: " << ms << "  FFInterp: " << FFInterp << "  Result: " << Result[i][f][n] << std::endl;
				}
				else
				{
					EmitResult[i][f][n] = 0;
				}
				
			}
		}
	}

	return EmitResult;
}



ScatteringRateMatrix LOPhononAbsScatRateCalc(FormFactorStruct LOPhononFF, PoissonResult PResult, KGridStruct KGrid, LOPhonStruct LOPhononParam, double TL)
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
	double DThetha = 2 * Pi / (NumTh - 1);

	//independent variable used for integration from0 to 2*Pi
	double Theta;

	//Momentum of LO Phonon involved in scattering
	double q;

	//Final Kf^2 of scattered electron
	double Kf;

	//Find Maximum energy of electron before leaving the well, Emax
	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());

	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Compute Scattering Rate Integral for initial subband i, final subband f, and magnitude of initial electron momentum ki

	//Initialize Phonon Emission Scattering Rate Matrix zero
	ScatteringRateMatrix AbsResult(NumSb, std::vector<std::vector<double>>(NumSb, std::vector<double>(Numk, 0)));

	//Initialize Ei and ms, Initial Energy and Effective mass
	double Ei, ms;

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
				//Compute initial guess for Effective Mass
				ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions[i], PResult.EigenEnergies[i]);

				//Initial energy of electron for non-parabolicity 
				Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms * ec);

				//Calculate Effective mass to account for non-parabolicity for given Ei
				ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions[i], Ei);

				//Find initial energy of electron accounting for non-parabolicity, i.e. using the proper ms 
				Ei = PResult.EigenEnergies[i] + hbar * hbar * KGrid.KMagVec[n] * KGrid.KMagVec[n] / (2 * ms * ec);

				if ( (PResult.EigenEnergies[f] - Ei <= LOGaAs.ELO) & (Ei <= Emax))
				{
					for (int th = 0; th < NumTh; th++)
					{

						//Theta the independent integration variable
						Theta = th * DThetha;

						//Find Kf the magnitude of the final K vector using: sqrt( ki^2 + (2ms)/h^2 * (Ei - Ef + ELO) )
						Kf = sqrt(KGrid.KMagVec[n] * KGrid.KMagVec[n] + 2 * ms * ec / (hbar * hbar) * (PResult.EigenEnergies[i] - PResult.EigenEnergies[f] + LOGaAs.ELO));

						//Find the magnitude of q the phonon momentum using: sqrt( ki^2 + kf^2 - 2*ki*kf*cos(Theta) )
						q = sqrt(KGrid.KMagVec[n] * KGrid.KMagVec[n] + Kf * Kf - 2 * KGrid.KMagVec[n] * Kf * cos(Theta));

						//Using the q calculated above, interpolate the FormFactor for LO Phonon Scattering
						FFInterp = FormFactorInterpolate(LOPhononFF, i, f, q);

						//Accumulate the Result of: ms*e^2*ELO/(8*pi*hbar^3)(1/Epinf - 1/Eps)*(NLO +1) integral( A(q)/q dTHeta from 0 to 2*pi)
						AbsResult[i][f][n] += ms * ec * ec * LOGaAs.ELO * ec / (8 * Pi * hbar * hbar * hbar * Ep0) * (1 / LOGaAs.Epinf - 1 / LOGaAs.Eps) * LOGaAs.NLO * FFInterp / q * DThetha;

					}
					//Printing Vals for Debug
					//std::cout << "NLO:  " << LOGaAs.NLO << std::endl;
					//std::cout << "i: " << i << "  f: " << f << "  Kf: " << Kf << "  Ki: " << KGrid.KMagVec[n] << "  q: " << q << "  Ei: " << PResult.EigenEnergies[i] << "  Ef: " << PResult.EigenEnergies[f] <<  std::endl;
				}
				else
				{
					AbsResult[i][f][n] = 0;
				}

			}
		}
	}

	return AbsResult;
}


ScatteringRateEEMatrix EEScatRateCalc(FormFactorEEStruct EEFF, PoissonResult PResult, KGridStruct KGrid, double TL, int Numq)
{
	//Find Maximum energy of electron before leaving the well, Emax
	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());

	//SubE, subband energies same as EigenEnergies, variable used for readibility
	std::vector<double> SubE = PResult.EigenEnergies;

	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Calculate Effective mass at Emax, as an upper bound
	double ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions.back(), Emax);

	//KiMax is found using parabolic band
	double Kimax = sqrt(Emax * 2 * ms * ec) / hbar;

	//Number of Subbands, used for convenience
	int NumSb = PResult.NewWaveFunctions.size();

	//Number of initial k-values to iterate over
	int Numk = KGrid.KMagVec.size();

	//Number of Theta values to integrate over
	int NumTh = Numq;

	//Number of Alpha values to integrate over
	int NumAlpha = Numq;

	//Delta Theta 
	double DThetha = 2 * Pi / (NumTh - 1);

	//Delta Alpha
	double DAlpha = 2 * Pi / (NumAlpha - 1);

	//Delta Kg
	double DKg=KGrid.KMagVec[1]- KGrid.KMagVec[0];

	//independent variable used for integration from0 to 2*Pi
	double Theta;

	//independent variable used for integration from0 to 2*Pi
	double Alpha;

	//Momentum exchange between initial and final states of primary electron
	double qxy;

	//Use a simple screening model for e-e scattering, where qxy = qxy + qs, qscreen is dependent on carrier denisty for a given subband and effective mass of subband 
	double qscreen;

	//q step, resolution in qvals used in the Form Factor Calculation
	double qs = EEFF.qvals[1] - EEFF.qvals[0];

	//nq , momentum index for the EE Form Factor Calculation
	int nq;

	//Prefactor for the Scattering Rate: e^4/(4*pi*hbar^3*(4*pi*Ep)^2)
	double PreFac = ec*ec*ec*ec/(64*Pi*Pi*Pi*hbar*hbar*hbar * Ep0 * Ep0 * PResult.NewZStruct.ZPermitivity[0] * PResult.NewZStruct.ZPermitivity[0]);


	//Compute Scattering Rate Integral for initial subband 1st electron i, final subband 1st electron f, initial subband 2nd electron g, 
	//  2nd electron final subband H, theta angle between difference of mometna of 1st and 2nd electron, Alpha angle between initial 1st and 2nd electron, 
	// kg the momentum of the 2nd electron

	//Initialize Electron Electron Rate Matrix zero
	//ScatteringRateMatrix EEResult(NumSb, std::vector<std::vector<double>>(NumSb, std::vector<double>(Numk, 0)));
	ScatteringRateEEMatrix EEResult(PResult.NewWaveFunctions.size(), std::vector<std::vector<std::vector<std::vector<double>>>>(PResult.NewWaveFunctions.size(), std::vector<std::vector<std::vector<double>>>(PResult.NewWaveFunctions.size(), std::vector<std::vector<double>>(PResult.NewWaveFunctions.size(), std::vector<double>(Numk, 0)))));

	//Initialize Ei Initial Energy, DKo2 k-vector corresponding to difference in initial and final energy of the electrons
	double Ei, DKo2;

	//Kig2 : Kig^2 , Kig is the diffence in intial momentum between electron 1 and 2
	double Kig2;

	//Kdelta : Small value that is Kig2/100, used because sqrt function yields small error that makes qxy imaginary for DKo2=0
	double Kdelta;

	//Probability of 2nd electron state being occupied, currently uses simple fermi-distribution
	double PKg;

	//Initialize Interpolated FormFactor Value
	double FFInterp;

	//Initialize Accumulated Scattering Rate over subbands g and H
	double SRAccum=0;

	//Debug variable
	double Val = 0;

	//Indexing over initial subband i
	for (int i = 0; i < NumSb; i++)
	{
		//Compute initial guess for Effective Mass
		ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions[i], PResult.EigenEnergies[i]);

		//Compute the contribution to screening, qscreen
		qscreen = ec*ec*ms/(2 * Ep0 * PResult.NewZStruct.ZPermitivity[0] * Pi *hbar*hbar) * FermiDist(SubE[i], PResult.u, TL);

		//Indexing over final subband f
		for (int f = 0; f < NumSb; f++)
		{
			//Indexing over magnitude of ki	
			for (int n = 0; n < Numk; n++)
			{
				//Indexing over initial subband g for 2nd electron
				for (int g = 0; g < NumSb; g++)
				{
					//Indexing over final subband H for 2nd electron
					for (int H = 0; H < NumSb; H++)
					{
						
						//For Debug						
						//std::cout << "i: " << i << "  f: " << f << "  n: " << n << "  g: " << g << "  H: " << H << std::endl;
					
						//k-vector corresponding to difference in initial and final energy of the electrons
						DKo2 = 4 * ms / (hbar * hbar) * ec * (SubE[i] + SubE[g] - SubE[f] - SubE[H]);

						//Indexing over magnitude of kg	
						for (int v = 0; v < Numk; v++)
						{
							//Calculate Probablity that state |g, kg> is occupied using simple fermi distribution, pass energy of 2nd electon, u-fermi level and Temp
							PKg = FermiDist(SubE[g] + hbar * hbar * KGrid.KMagVec[v] * KGrid.KMagVec[v] / (2 * ms * ec), PResult.u, TL);

							for (int a = 0; a < NumAlpha; a++)
							{
								//Theta the independent integration variable
								Alpha = a * DAlpha;

								//Calculate Kig^2 = Ki^2+Kg^2 -2*Ki*Kg*cos(alpha) 
								Kig2 = KGrid.KMagVec[n] * KGrid.KMagVec[n] + KGrid.KMagVec[v] * KGrid.KMagVec[v] - 2 * KGrid.KMagVec[n] * KGrid.KMagVec[v] * cos(Alpha);

								//Small value of k that prevents the calculation of qxy from being imaginary, due to round off error in sqrt
								Kdelta = sqrt(Kig2) / 1000;

								if ((DKo2 + Kig2 <= 0) or (Kig2 == 0))
								{
									// Integrate Nothing, Energy Conservation not conserved 
								}
								else
								{ 
									for (int th = 0; th < NumTh; th++)
									{
										//Theta the independent integration variable
										Theta = th * DThetha;

										//Find the magnitude of q the the momentum difference between Ki and Kf using: sqrt(2*Kig^2 + DKo^2 - 2*Kig*sqrt(Kig^2+DKo^2)*cos(Theta) ) /2, screen e-e interaction by adding qscreen
										qxy = sqrt(2 * Kig2 + DKo2 - 2 * sqrt(Kig2) * sqrt(Kig2 + DKo2) * cos(Theta) + Kdelta) / 2 + qscreen;

										//Using the qxy calculated above, interpolate the FormFactor for EE Scattering
										
										//Use index nq for the value of qxy, saves time compared to linear interpolation 
										nq = int(floor(qxy / qs));

										/*
										//For Debug purposes 
										if (n == 50)
										{
											std::cout << "i: " << i << "  f: " << f << "  g: " << g << "  H: " << H << " Kig2: " << Kig2 << "qxy: " << qxy << "  nq: " << nq << " theta: " << Theta << " Dko2: " << DKo2 << std::endl; // "   Upsampled FF: " << EEFF.FormFactor[i][f][g][H][nq] << std::endl;
										}
										
										//std::cout << "i: " << i << "  f: " << f << "  g: " << g << "  H: " << H << " Kig2: " << Kig2 << "qxy: " << qxy << "  nq: " << nq << std::endl; // "   Upsampled FF: " << EEFF.FormFactor[i][f][g][H][nq] << std::endl;
										*/
										
										//Accumulate the Result:
										//SRAccum += PKg * EEFF.FormFactor[i][f][g][H][nq] * EEFF.FormFactor[i][f][g][H][nq] /(qxy * qxy) * DThetha * DAlpha * KGrid.KMagVec[v] * DKg;
										
										Val = PKg * EEFF.FormFactor[i][f][g][H][nq] * EEFF.FormFactor[i][f][g][H][nq] / (qxy * qxy) * DThetha * DAlpha * KGrid.KMagVec[v] * DKg;

										SRAccum += Val;
										//std::cout << "i: " << i << "  f: " << f << "  g: " << g << "  H: " << H << "  Theta: " << Theta << "  Kn: " << KGrid.KMagVec[n] << "  Kv: " << KGrid.KMagVec[v] << "  Kig2: " << Kig2 << "  Dko2: " << DKo2 << " qxy: " << qxy << "  FF: " << EEFF.FormFactor[i][f][g][H][nq] << "  Val: " << Val << "  SRAccum: " << SRAccum  << std::endl; // "   Upsampled FF: " << EEFF.FormFactor[i][f][g][H][nq] << std::endl;
										

									}	
								}

							}
							
							//Printing Vals for Debug
							//std::cout << "NLO:  " << LOGaAs.NLO << std::endl;
							//std::cout << "i: " << i << "  f: " << f << "  Kf: " << Kf << "  Ki: " << KGrid.KMagVec[n] << "  q: " << q << "  Ei: " << PResult.EigenEnergies[i] << "  Ef: " << PResult.EigenEnergies[f] <<  std::endl;
						}
						// Accumulate the Scattering rate over subbands g and H, multiplied my the prefactor which includes effective mass
						EEResult[i][f][g][H][n] = PreFac * ms * SRAccum;

						//Reset SRAccum 
						SRAccum = 0;						

						//std::cout << "i: " << i << "  f: " << f << "  g: " << g << "  H: " << H << "	ki index  " << n  << "  Theta  " << Theta << "  Kig2: " << Kig2 << "  Dko2: " << DKo2 << " qxy: " << qxy << "  FF: " << EEFF.FormFactor[i][f][g][H][nq] << "  Val: " << Val << "  EERate: " << EEResult[i][f][g][H][n] << std::endl; // "   Upsampled FF: " << EEFF.FormFactor[i][f][g][H][nq] << std::endl;
					
						
						//std::cout << "EEResult:  " << EEResult[i][f][g][H][n] << std::endl;
						

					}
				}
				
			}
		} 
	}

	return EEResult;
}

