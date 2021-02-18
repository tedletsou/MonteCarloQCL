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


// Function to multiply vector by scalar
std::vector<double> MultiplyVectorByScalar(std::vector<double>& v, double k) 
{

    std::transform(v.begin(), v.end(), v.begin(), [k](double& c) { return c * k; });
    return v;
}
