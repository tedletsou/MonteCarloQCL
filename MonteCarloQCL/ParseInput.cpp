#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main() {

    // String for looping through .dat file
    string line;

    // Strings to compare .dat file values to
    string str1("laythick");
    string str2("laydop");
    string str3("laytype");
    string str4("modthick");
    string str5("field_vals");
    string str6("meshden");
    string str7("a_lat");
    string str8("mstar_bar");
    string str9("mstar_well");
    string str10("cband_bar");
    string str11("cband_well");
    string str12("vband_bar");
    string str13("vband_well");
    string str14("lhole_bar");
    string str15("lhole_well");
    string str16("splof_bar");
    string str17("splof_well");
    string str18("delso_bar");
    string str19("delso_well");
    string str20("Ep_bar");
    string str21("Ep_well");
    string str22("Eg_bar");
    string str23("Eg_well");

    // Inializing data variables

    vector<double> laythick;
    vector<double> laydop;
    vector<double> laytype;
    double modthick;
    vector<double> field_vals;
    double meshden;
    double a_lat;
    vector<double> mstar;
    vector<double> cband;
    vector<double> vband;
    vector<double> lhole;
    vector<double> sploff;
    vector<double> delso;
    vector<double> Ep;
    vector<double> Eg;

    // Load input file
    ifstream datafile("mcpp_input.dat");

    if (!datafile) {
        cout << "Unable to open file datafile.dat";
        exit(1);
    }

    while (datafile >> line) {
        // Layer thicknesses (monolayers)
        if (line.compare(str1) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            laythick.push_back(val);
        }
        // Layer dopings (cm^(-3))
        else if (line.compare(str2) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            laydop.push_back(val);
        }
        // Layer types (1: barrier, 2: well)
        else if (line.compare(str3) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            laytype.push_back(val);
        }
        // Module thickness (monolayers)
        else if (line.compare(str4) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            modthick = val;
        }
        // Applied fields (V/m)
        else if (line.compare(str5) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            field_vals.push_back(val);
        }
        // Meshing density
        else if (line.compare(str6) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            meshden = val;
        }
        // Lattice constant (monolayers)
        else if (line.compare(str7) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            a_lat = val;
        }
        // Effective mass in barrier and well (rel. to m0)
        else if (line.compare(str8) == 0 || line.compare(str9) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            mstar.push_back(val);
        }
        // Conduction band in barrier and well (eV)
        else if (line.compare(str10) == 0 || line.compare(str11) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            cband.push_back(val);
        }
        // Valence band in barrier and well (eV)
        else if (line.compare(str12) == 0 || line.compare(str13) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            vband.push_back(val);
        }
        // Light hole band in barrier and well (eV)
        else if (line.compare(str14) == 0 || line.compare(str15) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            lhole.push_back(val);
        }
        // Split-off band in barrier and well (eV)
        else if (line.compare(str16) == 0 || line.compare(str17) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            sploff.push_back(val);
        }
        // Delta split-off in barrier and well (eV)
        else if (line.compare(str18) == 0 || line.compare(str19) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            delso.push_back(val);
        }
        // Kane energies in barrier and well (eV)
        else if (line.compare(str20) == 0 || line.compare(str21) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            Ep.push_back(val);
        }
        // Band gap energies in barrier and well (eV)
        else if (line.compare(str22) == 0 || line.compare(str23) == 0) {
            getline(datafile, line, '\n');
            float val = stof(line);
            Eg.push_back(val);
        }
    }

    // Close input data file
    datafile.close();

    return 0;
}
