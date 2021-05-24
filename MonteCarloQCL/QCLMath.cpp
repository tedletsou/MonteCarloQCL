#include <vector>
#include <algorithm>
#include "QCLMath.h"
#include "Constants.h"
double FermiDist(double E, double u, double TL)
{
    double F = 1 / (1 + exp((E - u) * ec / (kb * TL)));

    return (F);
}

std::vector<double> Sinc1D(double Omega, int Length)
{
    std::vector<double> Y(Length, 0);

    int Midpoint = int((Length-1) / 2);

    for (int n = 0; n < Midpoint; n++)
    {
        Y[n] = sin(Omega * (n - Midpoint)) / (Omega * (n - Midpoint));
    }

    Y[Midpoint] = 1;

    for (int n = Midpoint+1; n < Length; n++)
    {
        Y[n] = sin(Omega * (n - Midpoint)) / (Omega * (n - Midpoint));
    }

    return(Y);
}


std::vector<double> Conv1D(std::vector<double> X, std::vector<double> H)
{
    int Nconv = int(H.size() + X.size() - 1);

    int HStart, XStart, XEnd;

    std::vector<double> Y(Nconv, 0);

    for (int i = 0; i < Nconv; i++)
    {
        XStart = std::max(0, int(i - H.size() + 1));
        XEnd = std::min(i + 1, int(X.size()));

        HStart = std::min(i, int(H.size() - 1));

        for (int j = XStart; j < XEnd; j++)
        {
            Y[i] += H[HStart--] * X[j];
        }
    }

    return Y;
}


// Linspace function similar to matlab
std::vector<double> linspace(double start_in, double end_in, int num_in)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
                              // are exactly the same as the input
    return linspaced;
}

//Linear Interpolation for 1D vector Upsampling, linear interpolation for upsampling
std::vector<double> UpSampleLin(std::vector<double> X, double UpFactor)
{
    //Initialize Output Vector Y to X
    std::vector<double> Y;

    std::vector<double> Temp(UpFactor, 0);

    //Linear Interpolate between values in X vector 
    for (int n = 0; n < X.size(); n++)
    {
        //Linear Inerpolation between adjacent points of X
        Temp = linspace(X[n], X[n + 1], UpFactor + 1);

        //Insert the Interpolated points except the last value into the upsampled vector Y
        Y.insert(Y.end(), Temp.begin(), Temp.end() - 1);
    }

    return(Y);
}




// Dot Product of Two Vectors
double operator*(const std::vector<double>& a, const std::vector<double>& x) {
    int m = a.size();

    double prod = 0;

    for (int i = 0; i < m; i++) {
        prod += a[i] * x[i];
    }
    return prod;
}

// Vector Multiply by a scalar
std::vector<double> operator*(double& a, const std::vector<double>& Vec)
{
    int m = Vec.size();

    std::vector<double> prod = Vec;

    for (int i = 0; i < m; i++) {
        prod[i] = a * Vec[i];
    }
    return prod;
}


// Vector Multiply by a scalar
std::vector<double> operator*(const std::vector<double>& Vec, double& a)
{
    int m = Vec.size();

    std::vector<double> prod = Vec;

    for (int i = 0; i < m; i++) {
        prod[i] = a * Vec[i];
    }
    return prod;
}

// Matrix Multiply by a Vector
std::vector<double> operator*(QCLMat& Mat, const std::vector<double>& Vec)
{
    int i, j;
    int m = Mat.size();
    int n = Vec.size();

    std::vector<double> prod(m);

    for (i = 0; i < m; i++) {
        prod[i] = 0.; 
        for (j = 0; j < n; j++)
            prod[i] += Mat[i][j] * Vec[j];
    }
    return prod;
}

// Matrix Multiply by a Matrix
QCLMat operator*(const QCLMat& a, const QCLMat& b) 
{
    int n = a.size();
    int m = a[0].size();
    int p = b[0].size();
    QCLMat c(n, std::vector<double>(p, 0));
    for (int j = 0; j < p; ++j) {
        for (int k = 0; k < m; ++k) {
            for (int i = 0; i < n; ++i) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    return c;
}

