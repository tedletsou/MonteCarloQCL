#include <vector>
#include <algorithm>
#include "QCLMath.h"

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

