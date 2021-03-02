#ifndef BOUNDEIGENENERGIES_H
#define BOUNDEIGENENERGIES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "GenerateZSpace.h"
#include "QCLMath.h"

QCLMat CalcEnergyBounds(ZMaterialParmsStruct ZStruct);

#endif