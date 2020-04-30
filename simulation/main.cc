#include "detectorrate.hh"
#include <sstream>
#include <iostream>
using namespace std;
int main(int argc, char *argv[]) {
    if (argc >= 4) {

        double Meff, Rbg, tdelay, binwidth;
        int fillopt;

        istringstream iss1(argv[1]);
        iss1 >> Meff; 
        istringstream iss2(argv[2]);
        iss2 >> Rbg;
        istringstream iss3(argv[3]);
        iss3 >> tdelay;

        if (argc >= 5) {
            istringstream iss4(argv[4]);
            iss4 >> binwidth;
        } else binwidth = 1;

        if (argc >= 6) {
            istringstream iss5(argv[5]);
            iss5 >> fillopt;
        } else fillopt = 1;

        detectorrate(Meff,Rbg,tdelay,binwidth,true,true,fillopt);

        return 0;
    }
    else {
        cerr << "run with: " << endl;
        cerr << argv[0] << " <Meff[kT]> <Rbg[Hz]> <tdelay[s]> (<bin width[ms] (<fillopt> 1 or 2) )" << endl;

    }
    return 1;
}
