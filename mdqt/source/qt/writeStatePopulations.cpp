#include "qt.h"

// this function saves the ion state populations vs time in a single file. the
// file format is given by [vx;statePop[0];...;statePop[end]];...repeat for each time point]
// , where each row has an entry for each ion and statePop[k] is the probability that
// a given ion occupies basis state |k>

using namespace std;

void QT::writeStatePopulations(std::filesystem::path savePath,const vector<double> &vx) const
// savePath: string containing save director for state populations relative to executable
// vx[i]: x-velocity for ion 'i' in units a*wPE
{
    // calculate p[k][i], the probability that ion 'i' occupies state 'k'
    vector<vector<double>> p(numStates);
    for (int k = 0; k < numStates; k++){
        p[k].resize(N0,0.);
    }
    for (int k = 0; k < numStates; k++){
        for (int i = 0; i < N0; i++){
            p[k][i] = pow(norm(bwvFn[k].t()*wvFn[i]),2.);
        }
    }

    // open state population file in append mode
    string delim{","};
    ofstream outFile(savePath,ofstream::app);

    // append ion velocities
    for (int i = 0; i < N0; i++){
        if (i < N0-1){
            outFile << vx[i] << delim;
        }
        else{
            outFile << vx[i] << endl;
        }
    }

    // append ion state populations
    for (int k = 0; k < numStates; k++){
        for (int i = 0; i < N0; i++){
            if (i < N0-1){
                outFile << p[k][i] << delim;
            }
            else{
                outFile << p[k][i] << endl;
            }
        }
    }

    // close state population file
    outFile.close();
}
