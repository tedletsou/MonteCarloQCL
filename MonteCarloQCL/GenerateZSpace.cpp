#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "QCLMath.h"
#include <fstream>

using namespace std;



int main() {

    // Temporary Vector used in loop to build Z Grid
    vector<double> TempVector;
    
    // Getting struct variables from parsed filed
    variables s = Parse("mcpp_input.dat");

    // naming outputs of struct for convenience
    ZThickness = s.laythick;
    MeshDen     = s.MeshDen;

    // sizing zspace to be used with partial_sum
    ZSpace.resize(ZThickness.size());

    // ZSpace is the cumulative thickness of each layer
    partial_sum(ZThickness.begin(), ZThickness.end(), ZSpace.begin());

    // Add a zero at beginning of ZSpace
    ZSpace.insert(ZSpace.begin(), 0);

    for (int k = 0; k < ZThickness.size(); k++)
    {
        //Temp Vector generates the grid for each layer
        TempVector = linspace(ZSpace[k], ZSpace[k + 1], MeshDen);

        // Remove the last element to remove redundant points at the layer boundary
        TempVector.pop_back();

        //Concatenate Zgrid with TempVector
        ZGrid.insert(ZGrid.end(),TempVector.begin(), TempVector.end());
    }

    
    for (int k = 0; k < ZThickness.size(); k++)
    {
        //Create Material Properties for wells and Barriers along Z, based on laytype (i.e. well or barrier material)
        CBE.insert(CBE.end(), MeshDen-1, s.cband[(int) s.laytype[k]-1 ]);
        ZMass.insert(ZMass.end(), MeshDen - 1, s.mstar[(int)s.laytype[k] - 1]);
        ZVBand.insert(ZVBand.end(), MeshDen - 1, s.vband[(int)s.laytype[k] - 1]);
        ZLHole.insert(ZLHole.end(), MeshDen - 1, s.lhole[(int)s.laytype[k] - 1]);
        ZSploff.insert(ZSploff.end(), MeshDen-1, s.sploff[(int)s.laytype[k] - 1]);
        ZDelta.insert(ZDelta.end(), MeshDen - 1, s.delso[(int)s.laytype[k] - 1]);
        ZKane.insert(ZKane.end(), MeshDen - 1, s.Ep[(int)s.laytype[k] - 1]);
        ZBandGap.insert(ZBandGap.end(), MeshDen - 1, s.Eg[(int)s.laytype[k] - 1]);
        
    }

    
    //Optional Code to write ZGrid and layer stucture to a File
    /*
    FILE* fpCBE = fopen("LayerStructure.txt", "w+");

    for (int k = 0; k < CBE.size(); k++)
    {
        fprintf(fpCBE, "%f \n", CBE[k]);

    }

    fclose(fpCBE);
    
    FILE* fp = fopen("ZGridCheck.txt", "w+");  

    for (int k = 0; k < ZGrid.size(); k++)
    {
        fprintf(fp, "%f \n", ZGrid[k]);

    }

    fclose(fp);
    */
    

    return 0;

}