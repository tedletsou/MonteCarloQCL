#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <complex>
#include <vector>
#include <stdio.h>
#include "EigenEnergyCalc.h"
#include "QCLMath.h"
#include "Constants.h"
#include "GenerateZSpace.h"
#include "EffectiveMass.h"


double ShootRoot(double WfEnergy, void* params)
{
    struct SolverParams* p = (struct SolverParams*)params;

    std::vector<double> Potential = p->Potential;
    ZMaterialParmsStruct ZStruct = p->ZStruct;


    //Following Physics of Photonic Devices 

    //Initial start of the wf, Psi(F)=0  Psi'(G)=1
    std::vector<double> F(Potential.size(), 0.0);
    std::vector<double> G(F.size(), 1.0);

    //Declare k vector
    std::complex<double> k(0, 0);

    //Calculate Effective Mass
    std::vector<double> ms = CalcEffectiveMass(WfEnergy, ZStruct);

    //Number of Nodes in the wavefunction, i.e. everytime F(z) changes sign
    int NumZero = 0;

    for (int n = 0; n < Potential.size() - 1; n++)
    {

        std::complex<double> ksquared = 2 * ms[n] * ec * (WfEnergy - Potential[n]) * pow(hbar, -2);

        //K-vector, sqrt(2ms(E-V)/hb^2)
        k = std::sqrt(ksquared);

        // Calculate difference in node position in meters
        double dz = (ZStruct.ZGridm[n + 1] - ZStruct.ZGridm[n]);

        //Calculate New Values of Psi (F) and Psi' (G) for next Z value

        // Taylor expand around k if value of k = 0;
        if (std::abs(k) != 0)
        {
            F[n + 1] = std::real(std::cos(k * dz) * F[n] + ((ms[n] / k) * std::sin(k * dz) * G[n]));
            G[n + 1] = std::real(-1.0 * (k / ms[n]) * std::sin(k * dz) * F[n] + (std::cos(k * dz) * G[n]));
        }
        else
        {
            F[n + 1] = std::real(std::cos(k * dz) * F[n] + ((ms[n] * dz) * G[n]));
            G[n + 1] = std::real(-1.0 * (k / ms[n]) * std::sin(k * dz) * F[n] + (std::cos(k * dz) * G[n]));
        }
    }

    double LastVal = F.back();

    return(LastVal);
}




std::vector<double> EigenEnergyCalc(QCLMat EigenEnergyBounds, ZMaterialParmsStruct ZStruct, std::vector<double> Potential, double EnergyTolerance)
{
    //Vector used to Store Eigen Energies
    std::vector<double> EigenEnergies{};

    //Struct used to contain parameters for the function ShootRoot
    SolverParams ShootRootParams{ Potential, ZStruct };

    //Updated bounds from FSolver used to calculate Solver Convergence
    double ELo;
    double EHi;

    //Temp Energy takes the intermediate value of Fsolve
    double TempEnergy;

    //Solver Status used to track when Solver has Converged
    int SolverStatus;

    //Pointers to fsolver object Solver and Function F
    const gsl_root_fsolver_type* T;
    gsl_root_fsolver* Solver;
    gsl_function F;

    //Use the Brent-Dekker method of fsolver in GNU Scientific Library (GSL) Library
    T = gsl_root_fsolver_brent;
    Solver = gsl_root_fsolver_alloc(T);

    F.function = &ShootRoot;
    F.params = &ShootRootParams;

    for (int k = 0; k < EigenEnergyBounds[0].size(); k++)
    {
        //Set Parameters for Fsolver with Energy Bounds for State k
        gsl_root_fsolver_set(Solver, &F, EigenEnergyBounds[0][k], EigenEnergyBounds[1][k]);

        do 
        {
            //Run FSolve and Check Convergence using gsl_root_test_interval
            SolverStatus = gsl_root_fsolver_iterate(Solver);
            TempEnergy=gsl_root_fsolver_root(Solver);
            ELo = gsl_root_fsolver_x_lower(Solver);
            EHi = gsl_root_fsolver_x_upper(Solver);
            SolverStatus = gsl_root_test_interval(ELo, EHi,0, EnergyTolerance);

            //Used to Debug convergence
            //printf(" [%.7f, %.7f] %.7f \n", ELo, EHi, TempEnergy);

        }
        while (SolverStatus == GSL_CONTINUE);
        EigenEnergies.push_back(TempEnergy);
    }

    gsl_root_fsolver_free(Solver);

    return(EigenEnergies);
}