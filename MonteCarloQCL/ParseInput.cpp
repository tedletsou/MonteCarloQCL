#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "ParseInput.h"

DeckDataStuct Parse(std::string fname) {

    // String for looping through .dat file
    std::string line;

    // Strings to compare .dat file values to
    std::string str1("laythick");
    std::string str2("laydop");
    std::string str3("laytype");
    std::string str4("modthick");
    std::string str5("field_vals");
    std::string str6("meshden");
    std::string str7("a_lat");
    std::string str8("mstar_bar");
    std::string str9("mstar_well");
    std::string str10("cband_bar");
    std::string str11("cband_well");
    std::string str12("vband_bar");
    std::string str13("vband_well");
    std::string str14("lhole_bar");
    std::string str15("lhole_well");
    std::string str16("sploff_bar");
    std::string str17("sploff_well");
    std::string str18("delso_bar");
    std::string str19("delso_well");
    std::string str20("Ep_bar");
    std::string str21("Ep_well");
    std::string str22("Eg_bar");
    std::string str23("Eg_well");
    std::string str24("Perm_bar");
    std::string str25("Perm_well");
    std::string str26("AlphaNP_bar");
    std::string str27("AlphaNP_well");
    std::string str28("BetaNP_bar");
    std::string str29("BetaNP_bar");

    // Load input file
    std::ifstream datafile(fname);

    if (!datafile) 
    {
        std::cout << "Unable to open file datafile.dat";
        exit(1);
    }

    DeckDataStuct s{};

    while (datafile >> line) {
        // Layer thicknesses (monolayers)
        if (line.compare(str1) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.laythick.push_back(val);
        }
        // Layer dopings (cm^(-3))
        else if (line.compare(str2) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.laydop.push_back(val);
        }
        // Layer types (1: barrier, 2: well)
        else if (line.compare(str3) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.laytype.push_back(val);
        }
        // Module thickness (monolayers)
        else if (line.compare(str4) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.modthick = val;
        }
        // Applied fields (V/m)
        else if (line.compare(str5) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.field_vals.push_back(val);
        }
        // Meshing density
        else if (line.compare(str6) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.MeshDen = val;
        }
        // Lattice constant (monolayers)
        else if (line.compare(str7) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.a_lat = val;
        }
        // Effective mass in barrier and well (rel. to m0)
        else if (line.compare(str8) == 0 || line.compare(str9) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.mstar.push_back(val);
        }
        // Conduction band in barrier and well (eV)
        else if (line.compare(str10) == 0 || line.compare(str11) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.cband.push_back(val);
        }
        // Valence band in barrier and well (eV)
        else if (line.compare(str12) == 0 || line.compare(str13) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.vband.push_back(val);
        }
        // Light hole band in barrier and well (eV)
        else if (line.compare(str14) == 0 || line.compare(str15) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.lhole.push_back(val);
        }
        // Split-off band in barrier and well (eV)
        else if (line.compare(str16) == 0 || line.compare(str17) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.sploff.push_back(val);
        }
        // Delta split-off in barrier and well (eV)
        else if (line.compare(str18) == 0 || line.compare(str19) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.delso.push_back(val);
        }
        // Kane energies in barrier and well (eV)
        else if (line.compare(str20) == 0 || line.compare(str21) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.Ep.push_back(val);
        }
        // Band gap energies in barrier and well (eV)
        else if (line.compare(str22) == 0 || line.compare(str23) == 0) 
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.Eg.push_back(val);
        }
        // Relative Permativity of barriers and well (unitless)
        else if (line.compare(str24) == 0 || line.compare(str25) == 0)
        {
            getline(datafile, line, '\n');
            float val = stof(line);
            s.Permitvity.push_back(val);
        }
        // Alpha term fon Non-Parabolicity of barriers and well (1/eV)
        else if (line.compare(str26) == 0 || line.compare(str27) == 0)
        {
        getline(datafile, line, '\n');
        float val = stof(line);
        s.alphaNP.push_back(val);
        }
        // Beta term fon Non-Parabolicity of barriers and well (1/eV)
        else if (line.compare(str28) == 0 || line.compare(str29) == 0)
        {
        getline(datafile, line, '\n');
        float val = stof(line);
        s.betaNP.push_back(val);
        }
    }

    // Close input data file
    datafile.close();

    // Returns the struct
    return s;
}