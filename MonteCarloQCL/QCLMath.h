#ifndef QCLMATH_H
#define QCLMATH_H

#include <vector>

// Defining linspace function (similar to MATLAB)
std::vector<double> linspace(double start_in, double end_in, int num_in);

// Defining vector scaling function
std::vector<double> MultiplyVectorByScalar(std::vector<double>& v, double k);

#endif
