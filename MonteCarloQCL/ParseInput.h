#ifndef PARSEINPUT_H
#define PARSEINPUT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Inializing data variables in struct
struct DeckDataStuct {
    std::vector<double> laythick;
    std::vector<double> laydop;
    std::vector<double> laytype;
    double modthick;
    std::vector<double> field_vals;
    double MeshDen;
    double a_lat;
    std::vector<double> mstar;
    std::vector<double> cband;
    std::vector<double> vband;
    std::vector<double> lhole;
    std::vector<double> sploff;
    std::vector<double> delso;
    std::vector<double> Ep;
    std::vector<double> Eg;
    std::vector<double> Permitvity;
};


// Defining parsing function
DeckDataStuct Parse(std::string fname);

#endif