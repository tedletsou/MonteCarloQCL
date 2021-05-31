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

		//Calculate initial Fermilevel u, Carrier Density in each subband and total charge density along Z (rho), returns only rho 
		std::vector<double> rho = CalcInitCarrierDensity(IonizedDopantDensity, ZMaterialStruct, WaveFunctions, EigenEnergies, TL);
		double u = CalcFermiLevel(IonizedDopantDensity, ZMaterialStruct, WaveFunctions, EigenEnergies, TL);

		//Poisson Solver Iterates to find no Change in Charge Density
		// !!!!!!!!!!!!!!!!! Hard Coded Error Tolerance for the Potential Convergence
		double ErrorTol = 1e0;
		PoissonResult PResult = PoissonSolver(ZMaterialStruct, IonizedDopantDensity, rho, WaveFunctions, TL, ErrorTol, DeckInput.field_vals[AppEindexa]);
		
		for (int n = 0; n < 1; n++)
		{
			std::cout << "New Z STruct Size" << PResult.NewZStruct.CBand.size() << std::endl;
		}

		//Calculate LO Phonon Occupation number, retrun LO phonon energy, permitivities and occupation number
		LOPhonStruct LOPhononParam = LOPhonGaAsOccupancy(TL);

		//Number of points for the Phonon Momentum =, Numq
		double Numq = 11;

		//Calculate the Form Factors for each subband for LO Phonon Scattering
		FormFactorStruct LOFF = FormFactorLOPhononCalc(PResult, LOPhononParam,Numq);

		std::cout << std::endl;
		std::cout << "Form Factor Calc Start" << std::endl;

		//Calculate the Form Factors for each subband for EE Scattering
		FormFactorEEStruct EEFF = FormFactorEECalc(PResult, Numq);

		std::cout << std::endl;
		std::cout << "Form Factor UpSample" << std::endl;

		// Upsample the EE Scatttering Form Factors by factor of Upfactor
		FormFactorEEStruct UpSampledEEFF = FormFactorUpSample(EEFF, 10);
		
		//Create Grid in K-Space, Kx, Ky, Kmag, Max K-Value is determined from Emax the engergy differnce from bottom to top of well
		KGridStruct KGrid = CreateKSpaceGrid(21, PResult, DeckInput);

		//Calculate LO Phonon Emission Scattering Rate in terms of Ki, i and f, the initial k-vector magnitude Ki, the initial subband i, the final subband f  
		ScatteringRateMatrix LOEmitScatRate = LOPhononEmitScatRateCalc(LOFF, PResult, KGrid, LOPhononParam, TL);

		//Calculate LO Phonon Absorption Scattering Rate in terms of Ki, i and f, the initial k-vector magnitude Ki, the initial subband i, the final subband f
		ScatteringRateMatrix LOAbsScatRate = LOPhononAbsScatRateCalc(LOFF, PResult, KGrid, LOPhononParam, TL);

		std::cout << std::endl << "Start your Engines!!" << std::endl;

		//Calculate EE Scattering Rate in terms of Ki, i and f, the initial k-vector magnitude Ki, the initial subband i, the final subband f
		ScatteringRateEEMatrix EEScatRate = EEScatRateCalc(UpSampledEEFF, PResult, KGrid, TL, Numq);



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

		//Optional Code to write out LO Phonon Form Factor to file
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

		//Optional Code to write out EE Form Factor to file
		FILE* fpFFEE = fopen("FormFactorEEOriginal.txt", "w+");

		fprintf(fpFFEE, "\n");

		//Indexing over initial subband i
		for (int i = 0; i < EEFF.FormFactor.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < EEFF.FormFactor[0].size(); f++)
			{
				//Indexing over final subband g
				for (int g = 0; g < EEFF.FormFactor[0][0].size(); g++)
				{
					//Indexing over final subband H
					for (int H = 0; H < EEFF.FormFactor[0][0][0].size(); H++)
					{
						//Indexing over magnitude of ki	
						for (int n = 0; n < EEFF.FormFactor[0][0][0][0].size(); n++)
						{
							fprintf(fpFFEE, "%d \t", i);
							fprintf(fpFFEE, "%d \t", f);
							fprintf(fpFFEE, "%d \t", g);
							fprintf(fpFFEE, "%d \t", H);
							fprintf(fpFFEE, "%.6g \t", EEFF.qvals[n]);
							fprintf(fpFFEE, "%.6g \t", EEFF.FormFactor[i][f][g][H][n]);
							fprintf(fpFFEE, "\n");
						}
					}
				}
			}
		}
		fclose(fpFFEE);

		//Optional Code to write out UpSampled EE Form Factor to file
		FILE* fpFFEEUP = fopen("FormFactorEEUpSampled.txt", "w+");

		fprintf(fpFFEEUP, "\n");

		//Indexing over initial subband i
		for (int i = 0; i < UpSampledEEFF.FormFactor.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < UpSampledEEFF.FormFactor[0].size(); f++)
			{
				//Indexing over final subband g
				for (int g = 0; g < UpSampledEEFF.FormFactor[0][0].size(); g++)
				{
					//Indexing over final subband H
					for (int H = 0; H < UpSampledEEFF.FormFactor[0][0][0].size(); H++)
					{
						//Indexing over magnitude of ki	
						for (int n = 0; n < UpSampledEEFF.FormFactor[0][0][0][0].size(); n++)
						{
							fprintf(fpFFEEUP, "%d \t", i);
							fprintf(fpFFEEUP, "%d \t", f);
							fprintf(fpFFEEUP, "%d \t", g);
							fprintf(fpFFEEUP, "%d \t", H);
							fprintf(fpFFEEUP, "%.6g \t", UpSampledEEFF.qvals[n]);
							fprintf(fpFFEEUP, "%.6g \t", UpSampledEEFF.FormFactor[i][f][g][H][n]);
							fprintf(fpFFEEUP, "\n");
						}
					}
				}
			}
		}
		fclose(fpFFEEUP);

		//Optional Code to write out LO Phonon Emission Rate to file
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
				
		//Optional Code to write out LO Phonon Absorption Rate to file
		FILE* fpSA = fopen("ScatteringRateLOAbsorb.txt", "w+");
		fprintf(fpSA, "\n");
		//Indexing over initial subband i
		for (int i = 0; i < LOAbsScatRate.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < LOAbsScatRate[0].size(); f++)
			{
				//Indexing over magnitude of ki	
				for (int n = 0; n < LOAbsScatRate[0][0].size(); n++)
				{
					fprintf(fpSA, "%d \t", i);
					fprintf(fpSA, "%d \t", f);
					fprintf(fpSA, "%.6g \t", KGrid.KMagVec[n]);
					fprintf(fpSA, "%.6g \t", LOAbsScatRate[i][f][n]);
					fprintf(fpSA, "\n");
				}

			}

		}
		fclose(fpSA);

		//Optional Code to write out EE Scattering Rate to file
		FILE* fpSEE = fopen("ScatteringRateEE.txt", "w+");
		fprintf(fpSEE, "\n");
		//Indexing over initial subband i
		for (int i = 0; i < EEScatRate.size(); i++)
		{
			//Indexing over final subband f
			for (int f = 0; f < EEScatRate[0].size(); f++)
			{
				//Indexing over final subband g
				for (int g = 0; g < EEScatRate[0][0].size(); g++)
				{
					//Indexing over final subband H
					for (int H = 0; H < EEScatRate[0][0][0].size(); H++)
					{
						//Indexing over magnitude of ki	
						for (int n = 0; n < EEScatRate[0][0][0][0].size(); n++)
						{
							fprintf(fpSEE, "%d \t", i);
							fprintf(fpSEE, "%d \t", f);
							fprintf(fpSEE, "%d \t", g);
							fprintf(fpSEE, "%d \t", H);
							fprintf(fpSEE, "%.6g \t", KGrid.KMagVec[n]);
							fprintf(fpSEE, "%.6g \t", EEScatRate[i][f][g][H][n]);
							fprintf(fpSEE, "\n");
						}

					}

				}
			}
		}
		fclose(fpSEE);
		
		
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