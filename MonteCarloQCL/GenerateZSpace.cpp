#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "ParseInput.h"

using namespace std;


int main() {

    variables s = Parse("mcpp_input.dat");
    cout << s.Eg[0] << endl;
    cout << s.Ep[0] << endl;
    return 0;

}