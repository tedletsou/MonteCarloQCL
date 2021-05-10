#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include "FormFactorCalc.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Shoot.h"
#include "PoissonSolver.h"
#include "EffectiveMass.h"
#include "Shoot.h"
#include "QCLMath.h"
#include "PhononPop.h"
#include "Constants.h"


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
	std::vector<double> msz = CalcEffectiveMass(Emax, PResult.NewZStruct);

	//Average over z, for the effective mass ms
	double ms = std::accumulate(msz.begin(), msz.end(), 0.0)/msz.size();

	//KiMax is found using parabolic band, sqrt(2) factor accounts for kmag^2=kx^2+ky^2
	double Kimax = sqrt(2) * sqrt(Emax * 2 * ms * ec) / hbar;

	double Kfmax = sqrt(Kimax * Kimax + 2 * ms * ec / (hbar * hbar) * (PResult.EigenEnergies.back() - PResult.EigenEnergies.front() + PhononParam.ELO));

	//qmin is 0
	double qmin = 0;

	//qmax is twice Kmax
	double qmax = Kimax + Kfmax;

	//For Debug Purpose

	std::cout << std::endl;
	std::cout << "Kimax value: " << Kimax << "Kfmax value: " << Kfmax << "Qmax value: " << qmax << std::endl;



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
