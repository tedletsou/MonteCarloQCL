#ifndef PARSEINPUT_H
#define PARSEINPUT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Inializing data variables in struct

struct variables {
    vector<double> laythick;
    vector<double> laydop;
    vector<double> laytype;
    double modthick;
    vector<double> field_vals;
    double MeshDen;
    double a_lat;
    vector<double> mstar;
    vector<double> cband;
    vector<double> vband;
    vector<double> lhole;
    vector<double> sploff;
    vector<double> delso;
    vector<double> Ep;
    vector<double> Eg;
};

variables Parse(string fname);

#endif