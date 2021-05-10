#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "BoundEigenDebug.h"
#include "ChargeDensityCalc.h"
#include "QCLMath.h"
#include "Shoot.h"
#include "BoundEigenenergies.h"
#include "EigenEnergyCalc.h"
#include "CalcWaveFunction.h"
#include "PoissonSolver.h"
#include "KGridSet.h"
#include "FormFactorCalc.h"
#include "PhononPop.h"
#include "ScattertingRateCalc.h" 

int main() 
{
	//Parse the Input Data for the simulation from mcpp_input.dat
	DeckDataStuct DeckInput = Parse("mcpp_input_thz.dat");

	//Create Vectors for the material parameters along the Z, Vectors stored in Struct ZMaterialParmsStruct
	ZMaterialParmsStruct ZMaterialStruct = CreateZParams(DeckInput);

	for (int AppEindexa =0; AppEindexa < DeckInput.field_vals.size(); AppEindexa++)
	{
		//Calculate the initial Potential of the QCL Structure with applied Bias and Conduction Band edge
		ZMaterialStruct = CalcPotential(ZMaterialStruct, DeckInput.field_vals[AppEindexa]);
		
		//!!!!!!!!! HARD CODED TEMP NEEDS to CHANGE
		double TL = 10;

		//Calculate initial Dopant Ion Distribution and Fermi Levels from Dopant Profile and Temperature 
		ChargeDistSturct IonizedDopantDensity = CalcInitDopantDensity(ZMaterialStruct, TL);

		// Find the Energy Bounds for each State in the Conduction Band
		QCLMat EnergyBounds = CalcEnergyBounds(ZMaterialStruct);

		
		//Allowed Error in Eigen Energies for Bound States
		double EnergyTolerance = 1e-8;

		// Find the Energy of each State in the Conduction Band based of Energy Bounds found above, using root finder in GNU Scientific Library (GSL), [Must inlcude in Project to function]
		std::vector<double> EigenEnergies = EigenEnergyCalc(EnergyBounds, ZMaterialStruct, EnergyTolerance);

		//Print out initial Eigen Energies
		std::cout << std::endl << "Initial Bound States  " << std::endl;
		for (int n = 0; n < EigenEnergies.size(); n++)
		{
			std::cout << EigenEnergies[n] << std::endl;
		}

		//Find Initial Wavefunctions based on Calculated Eigen Energies
		std::vector<WFStruct> WaveFunctions = CalculateWaveFunctions(EigenEnergies, ZMaterialStruct);

		//Calculate initial Fermilevel, Carrier Density in each subband and total charge density along Z (rho), returns only rho 
		std::vector<double> rho = CalcInitCarrierDensity(IonizedDopantDensity, ZMaterialStruct, WaveFunctions, EigenEnergies, TL);

		//Poisson Solver Iterates to find no Change in Charge Density
		// !!!!!!!!!!!!!!!!! Hard Coded Error Tolerance for the Potential Convergence
		double ErrorTol = 1e0;
		PoissonResult PResult = PoissonSolver(ZMaterialStruct, IonizedDopantDensity, rho, WaveFunctions, TL, ErrorTol, DeckInput.field_vals[AppEindexa]);
		
		for (int n = 0; n < 1; n++)
		{
			std::cout << "New Z STruct Size" << PResult.NewZStruct.CBand.size() << std::endl;
		}

		LOPhonStruct LOPhononParam = LOPhonGaAsOccupancy(TL);

		//Number of points for the Phonon Momentum =, Numq
		double Numq = 11;

		FormFactorStruct LOFF = FormFactorLOPhononCalc(PResult, LOPhononParam,Numq);

		KGridStruct KGrid = CreateKSpaceGrid(101, PResult, DeckInput);

		ScatteringRateMatrix LOEmitScatRate = LOPhononEmitScatRate(LOFF, PResult, KGrid, LOPhononParam, TL);


		/* Use to Bypass Poisson Solver and comment out Poissson Solver Call
		PoissonResult PResult;

		PResult.NewZStruct = ZMaterialStruct;
		PResult.NewWaveFunctions = WaveFunctions;
		*/

		//Optional Code to write EigenEnergies and Wavefunctions

		FILE* fpWF = fopen("WaveFunctions.txt", "w+");

		for (int n = 0; n < EigenEnergies.size(); n++)
		{
			fprintf(fpWF, "%.20g \t", EigenEnergies[n]);
		}

		fprintf(fpWF, "\n");

		for (int k = 0; k < WaveFunctions[0].Wavefunction.size(); k++)
		{
			for (int n = 0; n < WaveFunctions.size(); n++)
			{
				fprintf(fpWF, "%.20g \t", WaveFunctions[n].Wavefunction[k]);
			}
			fprintf(fpWF, "\n");
		}

		fclose(fpWF);


		//Optional Code to write Potential without Poisson Effect
		FILE* fpCBE = fopen("Potential.txt", "w+");

		for (int k = 0; k < ZMaterialStruct.Potential.size(); k++)
		{
			fprintf(fpCBE, "%f \t", ZMaterialStruct.Potential[k]);
			fprintf(fpCBE, "%f \n", ZMaterialStruct.ZGridm[k]*1e10);
		}

		fclose(fpCBE);

		//Optional Code to write Potential with Poisson Effect
		FILE* fpPB = fopen("PotentialBend.txt", "w+");

		for (int k = 0; k < PResult.NewZStruct.Potential.size(); k++)
		{
			fprintf(fpPB, "%f \t", PResult.NewZStruct.Potential[k]);
			fprintf(fpPB, "%f \n", PResult.NewZStruct.ZGridm[k] * 1e10);
		}

		fclose(fpPB);

		//Optional Code to write rho
		FILE* fpRoe = fopen("Rho.txt", "w+");

		for (int k = 0; k < rho.size(); k++)
		{
			fprintf(fpRoe, "%f \t", rho[k]);
			fprintf(fpRoe, "%f \n", PResult.NewZStruct.ZGridm[k] * 1e10);
		}

		fclose(fpRoe);

		FILE* fpSC = fopen("ScatteringRateLOEmission.txt", "w+");

		fprintf(fpSC, "\n");

		//Indexing over initial subband i
		for (int i = 0; i < LOEmitScatRate.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < LOEmitScatRate[0].size(); f++)
			{
				//Indexing over magnitude of ki	
				for (int n = 0; n < LOEmitScatRate[0][0].size(); n++)
				{
					fprintf(fpSC, "%d \t", i);
					fprintf(fpSC, "%d \t", f);
					fprintf(fpSC, "%.6g \t", KGrid.KMagVec[n]);
					fprintf(fpSC, "%.6g \t", LOEmitScatRate[i][f][n]);
					fprintf(fpSC, "\n");
				}

			}

		}

		fclose(fpSC);

		FILE* fpFF = fopen("FormFactorLOPhonon.txt", "w+");

		fprintf(fpFF, "\n");

		//Indexing over initial subband i
		for (int i = 0; i < LOFF.FormFactor.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < LOFF.FormFactor[0].size(); f++)
			{
				//Indexing over magnitude of ki	
				for (int n = 0; n < LOFF.FormFactor[0][0].size(); n++)
				{
					fprintf(fpFF, "%d \t", i);
					fprintf(fpFF, "%d \t", f);
					fprintf(fpFF, "%.6g \t", LOFF.qvals[n]);
					fprintf(fpFF, "%.6g \t", LOFF.FormFactor[i][f][n]);
					fprintf(fpFF, "\n");
				}

			}

		}

		fclose(fpFF);

		/*
		//Print out K Grid Values
		std::cout << std::endl << "KSpace Grid Kx  " << std::endl;
		for (int n = 0; n < KGrid.KGridKx.size(); n++)
		{
			for (int m = 0; m < KGrid.KGridKx.size(); m++)
			{
				std::cout << KGrid.KGridKMag[n][m] << std::endl;
				std::cout << std::endl << "Kx Row  " << std::endl;
			}
		}
		*/	
	};



return 0;
}