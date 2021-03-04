#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include "CalcWaveFunction.h"
#include "EffectiveMass.h"

double FermiLevel2DSheetDensity(double u, void* params)
{
    struct FermiSolverParams* p = (struct FermiSolverParams*)params;

    std::vector<double> EigenEnergies = p->EigenEnergies;
    double T = p->T;
    std::vector<double> ms = p->ms;
    double Nstot = p->Nstot;

    //Total Carrier Sheet Density 
    double Ns= 0.0;

    //Effective 2D carrier Density
    double Neff;

    for (int k = 0; k < EigenEnergies.size(); k++)
    {
        //Effective 2D carrier Density
        Neff = ms[k] * kb * T / (Pi * pow(hbar, 2));

        //Add the Contribution to the sheet density for subband indexed by k
        Ns += Neff * log(1.0 + exp((u - EigenEnergies[k]) / (kb / ec * T)));
    }

    return(Ns - Nstot);
}


ChargeDistSturct CalcInitDopantDensity(ZMaterialParmsStruct ZStruct, double TL)
{
    //!!!!!! USELESSS Function for now !!!!!!!

    // If you want to include Dopant Ionization, inlude in this function

    //Initialize Charge Dist and Fermi Level along Z
    ChargeDistSturct IonizedDopantDistZ{};

    IonizedDopantDistZ.FermiLevelZ.resize(ZStruct.ZGrid.size());
    IonizedDopantDistZ.ChargeDistZ.resize(ZStruct.ZGrid.size());

    //Calculate Fermi Level along Z for points in Z grid
    for (int k = 0; k < ZStruct.ZDoping.size(); k++)
    {
        if (ZStruct.ZDoping[k] != 0.0)
        {
            /* Used for calculating Donor Ionization. Commmented out instead simply assume all Donors are ionized
            //Effective Electron Density following Sze, Semiconductor Devices: Physics and Technology 2nd Edition, pg 35
            double Nc = 2 * pow(kb * TL * ZStruct.ZMass[k] / (2 * Pi * pow(hbar, 2)), 1.5);

            // Ef = Ec + Kb*T*ln( Nd / Nc), simple expression for Fermi level using Boltzmann Approx in the Extrinsic Region, violated for systems where carriers freeze out 
            IonizedDopantDistZ.FermiLevelZ[k] = ZStruct.CBand[k] + (kb/ec * TL) * std::log(ZStruct.ZDoping[k] / Nc);
            
            // n = Nc * Exp((Ef-Ec)/KbT), carrier density using Boltzmann Approx, violated for systems where carriers freeze out and at very high temp
            IonizedDopantDistZ.ChargeDistZ[k] = Nc* pow(e, ((IonizedDopantDistZ.FermiLevelZ[k] - ZStruct.CBand[k])) / (kb/ec * TL));
            std::cout << "FermiLevel  " << IonizedDopantDistZ.FermiLevelZ[k] << "  Nc   " << Nc << std::endl;
            */
            IonizedDopantDistZ.ChargeDistZ[k] = ZStruct.ZDoping[k];
            IonizedDopantDistZ.FermiLevelZ[k] = 0;
        }
        else 
        {
            // If not explicitly doped set the carrier density to 0, ignore thermal filling from valence band
            IonizedDopantDistZ.FermiLevelZ[k] = 0;
            IonizedDopantDistZ.ChargeDistZ[k] = 0;
        }

        
    }

    return IonizedDopantDistZ;
}


std::vector<double> CalcInitCarrierDensity(ChargeDistSturct IonizedDopantDensity, ZMaterialParmsStruct ZStruct, std::vector<WFStruct> WaveFunctions, std::vector<double> EigenEnergies, double TL)
{
    //Calculates the INITIAL electron distibution using a thermal distribution, SHOULD ONLY be used at the beginning of simulation before populations are found

    // HARD CODED VALUE for tolerance in finding the Fermi Level
    double FermiTolerance = 1e-5;

    // Fermi level
    double u;

    //Charge Density along Z
    std::vector<double> rho;

    //Resize rho to length of Z grid
    rho.resize(ZStruct.ZGrid.size());

    //Electron Density along Z
    std::vector<double> Nz;

    //Resize Nz to length of Z grid
    Nz.resize(ZStruct.ZGrid.size());

    //Electron Sheet density for each subband
    std::vector<double> Nes(EigenEnergies.size(),0);

    //Total Sheet Density for a single module equal to sum of Nes
    double Ns = 0;

    //Integrate the Dopant Density
    //DeltaZ found from ZStruct
    double DeltaZ;

    for (int k = 0; k < IonizedDopantDensity.ChargeDistZ.size()-1; k++)
    {
        DeltaZ = (ZStruct.ZGridm[k+1] - ZStruct.ZGridm[k]);
        Ns += DeltaZ * IonizedDopantDensity.ChargeDistZ[k];
    }

    std::cout << std::endl << "Sheet Density from Dopants: " << Ns << std::endl << std::endl;

    //Find bounds for the Fermi Level using the lowest subband energy as the min and half way between the lowest and highest subband engergy as the max
    double FermiLevelLow = -1; //EigenEnergies[0] - FermiTolerance * 10.0;
    double FermiLevelHigh = EigenEnergies.back() + FermiTolerance * 10.0 ;
    
    //Calculate Effective Mass for all EigenEnergies
    std::vector<double> ms = CalcWeightedEffectiveMass(ZStruct, WaveFunctions, EigenEnergies);

    //Struct used to contain parameters for the function ShootRoot
    FermiSolverParams FermiLevelParams{EigenEnergies, TL, ms, Ns};

    //Updated bounds from FSolver used to calculate Solver Convergence
    double ELow;
    double EHigh;

    //Calculate the Fermi Level using Fsolve from GSL library

    //Temp Energy takes the intermediate value of Fsolve
    double TemEnergy;

    //Solver Status used to track when Solver has Converged
    int SolverState;

    //Pointers to fsolver object Solver and Function F
    const gsl_root_fsolver_type* G;
    gsl_root_fsolver* Solverr;
    gsl_function H;


    //Use the Brent-Dekker method of fsolver in GNU Scientific Library (GSL) Library
    G = gsl_root_fsolver_brent;
    Solverr = gsl_root_fsolver_alloc(G);

    H.function = &FermiLevel2DSheetDensity;
    H.params = &FermiLevelParams;

    //Set Parameters for Fsolver with Energy Bounds for State k
    gsl_root_fsolver_set(Solverr, &H, FermiLevelLow, FermiLevelHigh);

    do
    {
    //Run FSolve and Check Convergence using gsl_root_test_interval
    SolverState = gsl_root_fsolver_iterate(Solverr);
    TemEnergy = gsl_root_fsolver_root(Solverr);
    ELow = gsl_root_fsolver_x_lower(Solverr);
    EHigh = gsl_root_fsolver_x_upper(Solverr);
    SolverState = gsl_root_test_interval(ELow, EHigh, 0, FermiTolerance);

    //Used to Debug convergence
    //printf(" [%.12f, %.12f] %.12f \n", ELow, EHigh, TemEnergy);

    } while (SolverState == GSL_CONTINUE);
   
    // Calculated fermi level
    u= TemEnergy;

    //Write out the fermi level
    std::cout << "Initial Fermi Level:    " << u << std::endl;

    //Free memory of fsolver
    gsl_root_fsolver_free(Solverr);

    //Effective 2D carrier Density
    double Neff;

    //Calculate the Sheet carrier density for each subband using fermi level
    for (int k = 0; k < EigenEnergies.size(); k++)
    {
        //Effective 2D carrier Density
        Neff = ms[k] * kb * TL / (Pi * pow(hbar, 2));

        //The sheet density for subband indexed by k
        Nes[k] = Neff * log(1.0 + exp((u - EigenEnergies[k]) / (kb / ec * TL)));
    }
    
    //NTot used to Check charge neutrality
    double NTot=0;

    //Used for spacing in integration
    double DZ;
   
    //Calculate the electron density using each subband Nes and the associated Wavefunction
    for (int m = 0; m < EigenEnergies.size(); m++)
    {
        for (int k = 0; k < ZStruct.ZGrid.size(); k++)
        {
            //Delta Z 
            DZ = ZStruct.ZGridm[k + 1] - ZStruct.ZGridm[k];

            //Add the contribution of each subband along z
            Nz[k] += Nes[m] * WaveFunctions[m].Wavefunction[k] * WaveFunctions[m].Wavefunction[k];

        }
        //Check charge neutrality NTot should be equal to Ns
        NTot += Nes[m];
    }

    //Ncheck is used to Check Charge neutrality, Ncheck should be equal to Ns and NTot
    double NCheck = 0;
    for (int k = 0; k < ZStruct.ZGrid.size(); k++)
    {
        DZ = ZStruct.ZGridm[k + 1] - ZStruct.ZGridm[k];
        NCheck += Nz[k] * DZ;
    }

    //Calculate rho the 3D charge density in z
    for (int k = 0; k < ZStruct.ZGrid.size(); k++)
    {
        rho[k]= ec*(IonizedDopantDensity.ChargeDistZ[k] - Nz[k]);
    }

    //Print values of Subband Carrier Density Nes
    std::cout << std::endl << "Subband Carrier Densities:  " << std::endl;
    for (int k = 0; k < Nes.size(); k++)
    {
        std::cout << k << "   " << Nes[k] << std::endl;
    }

    // Testing Only
    // Printing out Nes and NTot and Ncheck all values should be the same
    std::cout << std::endl << "Charge Neutrality Check   " << std::endl;
    printf("  [%.12f , %.12f, %.12f] \n", Ns, NTot, NCheck);
    std::cout << std::endl;

    //Check Charge Neutrality, sum of rho should be zero
    double rhosum = 0;
    for (int k = 0; k < rho.size(); k++)
    {
        DZ = ZStruct.ZGridm[k + 1] - ZStruct.ZGridm[k];
        rhosum += rho[k]*DZ;
    }

    std::cout << "Sum of rho   " << rhosum << std::endl << std::endl;

    return(rho);
}

