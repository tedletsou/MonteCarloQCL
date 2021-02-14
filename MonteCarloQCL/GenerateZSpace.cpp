#include <iostream>
#include <vector>
#include <numeric>
#include "ParseInput.h"
#include "GenerateZSpace.h"

using namespace std;



int main() {

    // Getting struct variables from parsed filed
    variables s = Parse("mcpp_input.dat");

    // naming outputs of struct for convenience
    thicknesses = s.laythick;
    meshden     = s.meshden;

    // sizing zspace to be used with partial_sum
    zspace.resize(thicknesses.size());

    partial_sum(thicknesses.begin(), thicknesses.end(), zspace.begin());

    return 0;

}