#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "GenerateZSpace.h"
#include "CalcWaveFunction.h"

ZMaterialParmsStruct CalcPotential(ZMaterialParmsStruct ZStruct, double AppliedField);

