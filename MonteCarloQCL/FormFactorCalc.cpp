#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "FormFactorCalc.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_vector.h>
#include "Shoot.h"
#include "PoissonSolver.h"
#include "EffectiveMass.h"
#include "Shoot.h"
#include "QCLMath.h"
#include "PhononPop.h"
#include "Constants.h"

FormFactorEEStruct FormFactorUpSample(FormFactorEEStruct FF, int UpFactor)
{
	//Number of Subbands used in Form Factor Calculation
	int NumSb = FF.FormFactor.size();

	//Number of q-values used in original(Down) Form Factor
	int NumqDown = FF.qvals.size();

	//Number of q-values used in UpSampled(Up) Form Factor
	int NumqUp = NumqDown * UpFactor;

	//Create q-vector for UpSampled(Up) Form Factor
	std::vector<double> qUp = linspace(FF.qvals.front(), FF.qvals.back(), NumqUp);

	//Initialize the UpSampled Form Factor 5D-Matrix to zero
	 FormFactorEEMatrix UpFormFactor(NumSb, std::vector<std::vector<std::vector<std::vector<double>>>>(NumSb, std::vector<std::vector<std::vector<double>>>(NumSb, std::vector<std::vector<double>>(NumSb, std::vector<double>(NumqUp, 0)))));


	//Iterate over the Form Factor, Seeding the Upsampled Form Factor with Original Values
	for (int i = 0; i < NumSb; i++)
	{
		for (int f = 0; f < NumSb; f++)
		{
			for (int g = 0; g < NumSb; g++)
			{
				for (int H = 0; H < NumSb; H++)
				{
					//Linear Interpolation across the q-vector of the Form Factor
					UpFormFactor[i][f][g][H] = UpSampleLin(FF.FormFactor[i][f][g][H], UpFactor);
				}
			}
		}
	}

	//Initialize Output Form Factor Struct
	FormFactorEEStruct UpFormFactorStruct;

	//Load Output Form Factor Struct
	UpFormFactorStruct.FormFactor = UpFormFactor;
	UpFormFactorStruct.qvals = qUp;

	return(UpFormFactorStruct);
}


double FormFactorInterpolate(FormFactorStruct FF, int i, int f, double q)
{
	
	//Implements a Monatonic Cubic Spline Interpolation of the Form Factor 

	//Pointer to Array for Form Factor Interpolation
	double* FFArray;
	double* QArray;

	//Size of Array is Number of qvals
	FFArray = new double[FF.qvals.size()];
	QArray = new double[FF.qvals.size()];

	//Copy Form Factor vector<double> to arrays
	std::copy(FF.FormFactor[i][f].begin(), FF.FormFactor[i][f].end(), FFArray);
	std::copy(FF.qvals.begin(), FF.qvals.end(), QArray);

	//For Debug Purpose
	/*
	std::cout << std::endl;
	std::cout << "qi value: " << q << std::endl;

	std::cout <<  std::endl;
	std::cout << "Qvals Array " << std::endl;

	for (int n = 0; n < FF.qvals.size(); n++)
	{
		std::cout << QArray[n] << std::endl;
	}
	std::cout << std::endl;

	std::cout << "FF Array " << std::endl;

	for (int n = 0; n < FF.qvals.size(); n++)
	{
		std::cout << FFArray[n] << std::endl;
	}
	std::cout << std::endl;
	*/

	//Declare gsl interpolate accelerator used to index and store data for gsl interpolation functions
	gsl_interp_accel* acc = gsl_interp_accel_alloc();

	//Declare gsl spline used for interpolating Form Factors 
	gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear, FF.qvals.size());

	//Calculate spline for interpolation
	gsl_spline_init(spline, QArray, FFArray, FF.qvals.size());

	//Compute the interpolated form factor using the initialized gsl_spline and the q value for the FormFactor
	double FFInterp = gsl_spline_eval(spline, q, acc);

	//Free Memory for gsl_spline objects and Arrays
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	delete[] FFArray;
	delete[] QArray;
	

	//Code Below implements Linear Interpolation, replaced with cubic spline in gsl library
	/*
	//Initialize Iterator used to find Index of Qmin and Qmax that bounds q
	std::vector<double>::iterator Itt;

	//Finding the Lower bound index of q in the vector of qvals, i.e. returning the iterator to qval right below q
	Itt = std::lower_bound(FF.qvals.begin(), FF.qvals.end(), q);

	//Initalize Qmin and Qmax that bounds q
	double Qmin, Qmax;

	//Set the value of QMax to the value of the Iterator pointer
	Qmax = *Itt;

	//Decrease the Iterator by one to Qmin the qval right below q
	std::advance(Itt, -1);
	Qmin = *Itt;

	//Find the Indicies of the vector qvals for Qmin and Qmax 
	int IndexMin = Itt - FF.qvals.begin();
	int IndexMax = IndexMin + 1;

	//Linear Interpolation between FF(Qmin) and FF(Qmax)
	double FFInterp = (q-Qmin)/(Qmax-Qmin)*FF.FormFactor[i][f][IndexMax]+(Qmax-q) / (Qmax - Qmin) * FF.FormFactor[i][f][IndexMin];
	*/

	return FFInterp;
}

double FormFactorEEInterpolate(FormFactorEEStruct FF, int i, int f, int g, int H, double qxy)
{

	//Implements a Monatonic Cubic Spline Interpolation of the Form Factor 

	//Pointer to Array for Form Factor Interpolation
	double* FFArray;
	double* QArray;

	//Size of Array is Number of qvals
	FFArray = new double[FF.qvals.size()];
	QArray = new double[FF.qvals.size()];

	//Copy Form Factor vector<double> to arrays
	std::copy(FF.FormFactor[i][f][g][H].begin(), FF.FormFactor[i][f][g][H].end(), FFArray);
	std::copy(FF.qvals.begin(), FF.qvals.end(), QArray);

	//Declare gsl interpolate accelerator used to index and store data for gsl interpolation functions
	gsl_interp_accel* acc = gsl_interp_accel_alloc();

	//Declare gsl spline used for interpolating Form Factors 
	gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear, FF.qvals.size());

	//Calculate spline for interpolation
	gsl_spline_init(spline, QArray, FFArray, FF.qvals.size());

	//Compute the interpolated form factor using the initialized gsl_spline and the q value for the FormFactor
	double FFInterp = gsl_spline_eval(spline, qxy, acc);

	//Free Memory for gsl_spline objects and Arrays
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	delete[] FFArray;
	delete[] QArray;

	return FFInterp;
}


FormFactorStruct FormFactorLOPhononCalc(PoissonResult PResult, LOPhonStruct PhononParam, int Numq)
{
	PResult.NewZStruct.ZMass;
	// Find Maximum Phonon momentum qmax, using the maximum momentum transfer between -Kmax to Kmax, 
	// where Kmax is found from Emax the bounds of the Potential, Emax=max(Potential)-min(Potential)
	
	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());
	
	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Calculate Effective mass at Emax, as an upper bound
	double ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions.back(), Emax);

	//KiMax is found using parabolic band
	double Kimax = sqrt(Emax * 2 * ms * ec) / hbar;

	double Kfmax = sqrt(Kimax * Kimax + 2 * ms * ec / (hbar * hbar) * (Emax + PhononParam.ELO));

	//qmin is 0
	double qmin = 0;

	//qmax is twice Kmax
	double qmax = Kimax + Kfmax;

	//For Debug Purpose
	//std::cout << std::endl;
	//std::cout << "Emax value: " << Emax  << "  Kimax value: " << Kimax << "  Kfmax value: " << Kfmax << "  Qmax value: " << qmax << "  ms value: " << ms << std::endl;

	//Dq is the incremental q from qmin to qmax, qmin:Dq:qmax
	double Dq = (qmax - qmin) / (Numq-1);

	//Initialize q vector using linspace
	std::vector<double> q = linspace(qmin, qmax, Numq);

	//Initialize the Form Factor 3D-Matrix to zero
	FormFactorMatrix FormFactorLO(PResult.NewWaveFunctions.size(), std::vector<std::vector<double>>(PResult.NewWaveFunctions.size(), std::vector<double>(Numq, 0)));

	//Number of Subbands, used for convenience
	int NumSb = PResult.NewWaveFunctions.size();

	//Integrate WavefunctionI(z) * WavefunctionF(z) * WavefunctionI(z') * WavefunctionF(z') * exp(q*|z-z'|)
	//Delta Z
	double DZ;
	//Delta Zp
	double DZp;

	WFStruct WaveFunctionSt = PResult.NewWaveFunctions[0];

	//Temp Result
	std::vector<double> IntResult(Numq, 0);

	std::vector<double> Z = PResult.NewZStruct.ZGridm;
	std::vector<double> Zp = PResult.NewZStruct.ZGridm;
	
	//Indexing over initial subband i
	for (int i = 0; i < NumSb; i++)
	{
		//Indexing over final subband f
		for (int f = i; f < NumSb; f++)
		{
			//Indexing over q	
			for (int n = 0; n < Numq; n++)
			{
				//Indexing over z
				for (int k = 0; k < PResult.NewZStruct.ZGrid.size(); k++)
				{
					//Delta Z used for integration
					DZ = Z[k + 1] - Z[k];

					//Indexing over zp
					for (int m = 0; m < PResult.NewZStruct.ZGrid.size(); m++)
					{
						//Delta Zp used for integration
						DZp = Zp[m + 1] - Zp[m];

						FormFactorLO[i][f][n] += PResult.NewWaveFunctions[i].Wavefunction[m] * PResult.NewWaveFunctions[f].Wavefunction[m] * PResult.NewWaveFunctions[i].Wavefunction[k] * PResult.NewWaveFunctions[f].Wavefunction[k] * exp(-1 * abs(Z[k] - Zp[m]) * q[n]) * DZ * DZp;
					}
				}
			}
		}
	}

	// Populating Matrix Symetric Elements not calculated above, Optimization to avoid redundant calculation
	//Indexing over final subband f
	for (int f = 0; f < NumSb; f++)
	{
		//Indexing over initial subband i
		for (int i = f; i < NumSb; i++)
		{
			FormFactorLO[i][f] = FormFactorLO[f][i];
		}
	}
	
	/*
	//Print Out for Debug Purpose
	std::cout << "Qmax = " << qmax << "  Dq: " << Dq<< std::endl;

	//Print Out FormFactor for Debug Purpose
	for (int i = 0; i < NumSb; i++)
	{
		for (int f = 0; f < NumSb; f++)
		{
			std::cout << "Init Subband: " << i << "   Final Subband: " << f << std::endl;
			for (int n = 0; n < Numq; n++)
			{
				std::cout << q[n] << "  :  " << FormFactorLO[i][f][n] << std::endl;
			}
			std::cout << std::endl;
		}
	}
	*/

	FormFactorStruct FormFactorLOStruct;

	FormFactorLOStruct.FormFactor = FormFactorLO;
	FormFactorLOStruct.qvals = q;


	return FormFactorLOStruct;
}

FormFactorEEStruct FormFactorEECalc(PoissonResult PResult, int Numq)
{
	PResult.NewZStruct.ZMass;
	// Find Maximum Phonon momentum qmax, using the maximum momentum transfer between -Kmax to Kmax, 
	// where Kmax is found from Emax the bounds of the Potential, Emax=max(Potential)-min(Potential)

	// Find min and max of conduction band energy after Poisson Effect 
	auto MinMax = std::minmax_element(PResult.NewZStruct.Potential.begin(), PResult.NewZStruct.Potential.end());

	//Maximum Bound State Energy
	double Emax = *MinMax.second - *MinMax.first;

	//Calculate Effective mass at Emax, as an upper bound
	double ms = CalcAverageEffectiveMass(PResult.NewZStruct, PResult.NewWaveFunctions.back(), Emax);

	//KiMax is found using parabolic band
	double Kimax = sqrt(Emax * 2 * ms * ec) / hbar;

	double Kfmax = sqrt(Kimax * Kimax + 2 * ms * ec / (hbar * hbar) * (Emax));

	//qmin is slightly above zero to avoid divide by 0
	double qmin = 1e5;

	//qmax is twice Kmax
	double qmax = Kimax + Kfmax + qmin;

	//For Debug Purpose
	//std::cout << std::endl;
	//std::cout << "Emax value: " << Emax  << "  Kimax value: " << Kimax << "  Kfmax value: " << Kfmax << "  Qmax value: " << qmax << "  ms value: " << ms << std::endl;

	//Dq is the incremental q from qmin to qmax, qmin:Dq:qmax
	double Dq = (qmax - qmin) / (Numq - 1);

	//Initialize q vector using linspace
	std::vector<double> q = linspace(qmin, qmax, Numq);

	//Initialize the Form Factor 5D-Matrix to zero
	FormFactorEEMatrix FormFactorEE(PResult.NewWaveFunctions.size(), std::vector<std::vector<std::vector<std::vector<double>>>>(PResult.NewWaveFunctions.size(), std::vector<std::vector<std::vector<double>>>(PResult.NewWaveFunctions.size(), std::vector<std::vector<double>>(PResult.NewWaveFunctions.size(), std::vector<double>(Numq,0)))));

	//Number of Subbands, used for convenience
	int NumSb = PResult.NewWaveFunctions.size();

	//Integrate WavefunctionI(z) * WavefunctionF(z) * WavefunctionG(z') * WavefunctionH(z') * exp(q*|z-z'|)
	//Delta Z
	double DZ;
	//Delta Zp
	double DZp;

	WFStruct WaveFunctionSt = PResult.NewWaveFunctions[0];

	//Temp Result
	std::vector<double> IntResult(Numq, 0);

	std::vector<double> Z = PResult.NewZStruct.ZGridm;
	std::vector<double> Zp = PResult.NewZStruct.ZGridm;

	//Indexing over initial subband i
	for (int i = 0; i < NumSb; i++)
	{
		//Indexing over final subband f
		for (int f = i; f < NumSb; f++)
		{
			//Indexing over initial subband g for 2nd electron
			for (int g = 0; g < NumSb; g++)
			{
				//Indexing over final subband H for 2nd electron
				for (int H = g; H < NumSb; H++)
				{
					//Indexing over q	
					for (int n = 0; n < Numq; n++)
					{
						//Indexing over z
						for (int k = 0; k < PResult.NewZStruct.ZGrid.size(); k++)
						{
							//Delta Z used for integration
							DZ = Z[k + 1] - Z[k];

							//Indexing over zp
							for (int m = 0; m < PResult.NewZStruct.ZGrid.size(); m++)
							{
								//Delta Zp used for integration
								DZp = Zp[m + 1] - Zp[m];

								FormFactorEE[i][f][g][H][n] += PResult.NewWaveFunctions[i].Wavefunction[m] * PResult.NewWaveFunctions[f].Wavefunction[m] * PResult.NewWaveFunctions[g].Wavefunction[k] * PResult.NewWaveFunctions[H].Wavefunction[k] * exp(-1 * abs(Z[k] - Zp[m]) * q[n]) * DZ * DZp;
							}
						}
						// Compute F(i,f,g,h,n,) = FF^2/qxy^2, speeds up calculation of EE scattering rate in ScatteringRateCalc.cpp
						//FormFactorEE[i][f][g][H][n] = FormFactorEE[i][f][g][H][n] * FormFactorEE[i][f][g][H][n] / (q[n] * q[n]);
					}
				}
			}
		}
	}

// Populating Matrix Symetric Elements not calculated above, Optimization to avoid redundant calculation
	//Indexing over final subband f
	for (int f = 0; f < NumSb; f++)
	{
		//Indexing over initial subband i
		for (int i = f; i < NumSb; i++)
		{
			//Indexing over initial subband H
			for (int H = 0; H < NumSb; H++)
			{
				//Indexing over final subband g
				for (int g = H; g < NumSb; g++)
				{
					FormFactorEE[i][f][g][H] = FormFactorEE[f][i][H][g];
				}
			}
		}
	}

	/*
	//Print Out for Debug Purpose
	std::cout << "Qmax = " << qmax << "  Dq: " << Dq<< std::endl;

	//Print Out FormFactor for Debug Purpose
	for (int i = 0; i < NumSb; i++)
	{
		for (int f = 0; f < NumSb; f++)
		{
			for (int g = 0; g < NumSb; g++)
			{
				for (int H = 0; H < NumSb; H++)
				{
					std::cout << "Init Subband: " << i << "   Final Subband: " << f << "   2nd Initial Subband: " << g << "   2nd Final Subband: " << H << std::endl;
					for (int n = 0; n < Numq; n++)
					{
						std::cout << q[n] << "  :  " << FormFactorEE[i][f][g][H][n] << std::endl;
					}
					std::cout << std::endl;
				}
			}
		}
	}
	*/

	FormFactorEEStruct FormFactorStruct;

	FormFactorStruct.FormFactor = FormFactorEE;
	FormFactorStruct.qvals = q;


	return FormFactorStruct;
}