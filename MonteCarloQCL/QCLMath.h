#ifndef QCLMATH_H
#define QCLMATH_H

#include <vector>

typedef std::vector<double> QCLVec;

// Defining linspace function (similar to MATLAB)
std::vector<double> linspace(double start_in, double end_in, int num_in);

typedef std::vector<double> QCLVec;
typedef std::vector<std::vector<double>> QCLMat;

//double operator*(const QCLVec& a, const QCLVec& x);

double operator*(const std::vector<double>& a, const std::vector<double>& x);

std::vector<double> operator*(double& a, const std::vector<double>& Vec);

std::vector<double> operator*(QCLMat& Mat, const std::vector<double>& Vec);

std::vector<double> operator*(const std::vector<double>& Vec, double& a);

QCLMat operator*(const QCLMat& a, const QCLMat& b);


#endif

