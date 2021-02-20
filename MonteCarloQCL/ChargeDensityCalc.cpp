#include <iostream>
#include <fstream>
#include "ChargeDensityCalc.h"
#include "Constants.h"
#include "QCLMath.h"
#include <vector>

ChargeDistSturct CalcInitDopantDensity(ZMaterialParmsStruct ZStruct, double TL)
{
    //Initialize Charge Dist and Fermi Level along Z
    ChargeDistSturct InitChargeDistZ{};

    InitChargeDistZ.FermiLevelZ.resize(ZStruct.ZGrid.size());
    InitChargeDistZ.ChargeDistZ.resize(ZStruct.ZGrid.size());

    //Calculate Fermi Level along Z for points in Z grid
    for (int k = 0; k < ZStruct.ZDoping.size(); k++)
    {
        if (ZStruct.ZDoping[k] != 0)
        {
            //Effective Electron Density following Sze, Semiconductor Devices: Physics and Technology 2nd Edition, pg 35
            double Nc = 2 * pow(kb * TL * ZStruct.ZMass[k] / (2 * Pi * pow(hbar, 2)), 1.5);

            // Ef = Ec + Kb*T*ln( Nd / Nc), simple expression for Fermi level using Boltzmann Approx in the Extrinsic Region, violated for systems where carriers freeze out 
            InitChargeDistZ.FermiLevelZ[k] = ZStruct.CBand[k] + (kb * TL) * std::log(ZStruct.ZDoping[k] / Nc);
            
            // n = Nc * Exp((Ef-Ec)/KbT), carrier density using Boltzmann Approx, violated for systems where carriers freeze out and at very high temp
            InitChargeDistZ.ChargeDistZ[k] = Nc* pow(e, ((InitChargeDistZ.FermiLevelZ[k] - ZStruct.CBand[k])) / (kb * TL));
        }
        else 
        {
            // If not explicitly doped set the carrier density to 0, ignore thermal filling from valence band
            InitChargeDistZ.FermiLevelZ[k] = 0;
        }
    }

    return InitChargeDistZ;
}
