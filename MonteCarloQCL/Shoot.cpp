#include <math.h>
#include "Shoot.h"
#include "QCLMath.h"
#include "Constants.h"
#include <complex> 


WFStruct Shoot(double WfEnergy, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms)
{
    //Following Physics of Photonic Devices 

    //Declare Output Struct
    WFStruct Result;

    //Initial start of the wf, Psi(F)=0  Psi'(G)=1
    std::vector<double> F(Potential.size(), 0.0);
    std::vector<double> G(F.size(), 1.0);
    
    //Declare k vector
    std::complex<double> k (0, 0);

    //Number of Nodes in the wavefunction, i.e. everytime F(z) changes sign
    int NumZero=0;

    for (int n = 0; n < Potential.size()-1; n++)
    {
        
        std::complex<double> ksquared = 2 * ms[n] * ec * (WfEnergy - Potential[n]) * pow(hbar, -2);

        //K-vector, sqrt(2ms(E-V)/hb^2)
        k = std::sqrt(ksquared);

        // Calculate difference in node position in meters
        double dz = (ZGridm[n + 1] - ZGridm[n]);

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

        //Check to see if sign of F changed by summing NumZero with sign(F[n])-sign(F[n+1])
        NumZero += int(signbit(F[n]) ^ signbit(F[n + 1]));
    }

    
    Result.Wavefunction = F;
    Result.NumZeros = NumZero;
    Result.LastValue = F.back();
       
    return(Result);
}
