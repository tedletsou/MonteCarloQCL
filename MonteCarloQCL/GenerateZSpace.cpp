#pragma warning(disable : 4996)

#include <iostream>
#include <vector>
#include <numeric>
#include "ParseInput.h"
#include "GenerateZSpace.h"
#include "Constants.h"
#include "QCLMath.h"
#include <fstream>
#include <functional> 
#include <algorithm>

ZMaterialParmsStruct CreateZParams(DeckDataStuct DeckData) 
{

    //Thickness of Each layer along Z
    std::vector<double> ZThickness;

    //Intereface Posistions along the Z axis
    std::vector<double> ZSpace;

    //Create Struct to Load Z Parameters into
    ZMaterialParmsStruct ZStruct{};

    // Temporary Vector used in loop to build Z Grid
    std::vector<double> TempVector;

    // Temporary doping vector used in loop to build Z Grid
    std::vector<double> TempDopingVector;

    // naming outputs of struct for convenience
    ZThickness = DeckData.laythick;
    //MeshDensity = DeckData.MeshDen;

    // sizing zspace to be used with partial_sum
    ZSpace.resize(ZThickness.size());

    // ZSpace is the cumulative thickness of each layer
    partial_sum(ZThickness.begin(), ZThickness.end(), ZSpace.begin());

    // Add a zero at beginning of ZSpace
    ZSpace.insert(ZSpace.begin(), 0);

    for (int k = 0; k < ZThickness.size(); k++)
    {
        //Temp Vector generates the grid for each layer
        TempVector = linspace(ZSpace[k], ZSpace[k + 1], DeckData.MeshDen);

        // Temp vector for doping at each layer
        TempDopingVector = linspace(DeckData.laydop[k], DeckData.laydop[k], DeckData.MeshDen - 1);

        // Remove the last element to remove redundant points at the layer boundary
        TempVector.pop_back();

        //Concatenate Zgrid with TempVector
        ZStruct.ZGrid.insert(ZStruct.ZGrid.end(), TempVector.begin(), TempVector.end());

        // Concatenate ZDoping with TempDopingVector
        ZStruct.ZDoping.insert(ZStruct.ZDoping.end(), TempDopingVector.begin(), TempDopingVector.end());
    }

    //Calculate Zgrid in meters, ZGrid * lattice constant/2
    double MetersConversion = DeckData.a_lat * 0.5 * 1E-10;

    ZStruct.ZGridm = ZStruct.ZGrid * MetersConversion;
    
    // Constant to convert from cm^(-3) to m^(-3)
    double DopingConversion = 1E6;

    // Converting doping array from cm^(-3) to m^(-3)
    ZStruct.ZDoping = ZStruct.ZDoping * DopingConversion;
    
    for (int k = 0; k < ZThickness.size(); k++)
    {
        //Create Material Properties for wells and Barriers along Z, based on laytype (i.e. well or barrier material)
        ZStruct.CBand.insert(ZStruct.CBand.end(), DeckData.MeshDen -1, DeckData.cband[(int)DeckData.laytype[k]-1 ]);
        ZStruct.ZMass.insert(ZStruct.ZMass.end(), DeckData.MeshDen - 1, DeckData.mstar[(int)DeckData.laytype[k] - 1] * me);
        ZStruct.ZVBand.insert(ZStruct.ZVBand.end(), DeckData.MeshDen - 1, DeckData.vband[(int)DeckData.laytype[k] - 1]);
        ZStruct.ZLHole.insert(ZStruct.ZLHole.end(), DeckData.MeshDen - 1, DeckData.lhole[(int)DeckData.laytype[k] - 1]);
        ZStruct.ZSploff.insert(ZStruct.ZSploff.end(), DeckData.MeshDen -1, DeckData.sploff[(int)DeckData.laytype[k] - 1]);
        ZStruct.ZDelta.insert(ZStruct.ZDelta.end(), DeckData.MeshDen - 1, DeckData.delso[(int)DeckData.laytype[k] - 1]);
        ZStruct.ZKane.insert(ZStruct.ZKane.end(), DeckData.MeshDen - 1, DeckData.Ep[(int)DeckData.laytype[k] - 1]);
        ZStruct.ZBandGap.insert(ZStruct.ZBandGap.end(), DeckData.MeshDen - 1, DeckData.Eg[(int)DeckData.laytype[k] - 1]);
        
    }

    //Optional Code to write ZGrid and layer stucture to a File
    /*
    FILE* fpCBE = fopen("LayerStructure.txt", "w+");

    for (int k = 0; k < ZStruct.CBand.size(); k++)
    {
        fprintf(fpCBE, "%f \n", CBE[k]);

    }

    fclose(fpCBE);
    
    FILE* fp = fopen("ZDopingCheck.txt", "w+");  

    for (int k = 0; k < ZStruct.ZDoping.size(); k++)
    {
        fprintf(fp, "%f \n", ZStruct.ZDoping[k]);

    }

    fclose(fp);

    */

    //cout << ZStruct.ZDoping[130];

    return ZStruct;

}