#include <math.h>
#include "Shoot.h"
#include "QCLMath.h"
#include "Constants.h"


WFStruct Shoot(double WfEnergy, std::vector<double> ZGridm, std::vector<double> Potential, std::vector<double> ms)
{
    //Following Physics of Photonic Devices 

    //Declare Output Struct
    WFStruct Result;

    //Initial start of the wf, Psi(F)=0  Psi'(G)=1
    std::vector<double> F(Potential.size(),0);
    std::vector<double> G(Potential.size(), 1);
    
    //Declare k vector
    double k;

    //Number of Nodes in the wavefunction, i.e. everytime F(z) changes sign
    int NumZero=0;

    for (int n = 0; n < Potential.size()-1; n++)
    {
        //K-vector, sqrt(2ms(E-V)/hb^2)
        k = sqrt(2 * ms[n] / pow(hbar, 2) * ec * (WfEnergy - Potential[n]));

        //Calculate New Values of Psi (F) and Psi' (G) for next Z value
        F[n + 1] = cos(k * ZGridm[n]) * F[n] + (ms[k] / k) * sin(k * ZGridm[n]) * G[n];
        G[n + 1] = cos(k * ZGridm[n]) * G[n] - (ms[k] / k) * sin(k * ZGridm[n]) * F[n];

        //Check to see if sign of F changed by summing NumZero with sign(F[n])-sign(F[n+1])
        NumZero += int(signbit(F[n]) ^ signbit(F[n + 1]));
    }

    
    Result.Wavefunction = F;
    Result.NumZeros = NumZero;
    Result.LastValue = F.back();
       
    return(Result);
}
