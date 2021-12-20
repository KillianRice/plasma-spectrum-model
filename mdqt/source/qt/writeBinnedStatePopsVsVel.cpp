#include "qt.h"

void QT::writeBinnedStatePopsVsVel(std::filesystem::path filePath,vector<double> vBin, vector<double> vx,double time) const
// filePath: output file save path relative to runFile
// vBin[j]: x-velocity (in SI units) for velocity bin 'j'
// vx[i]: x-velocity for ion 'i' (in SI units)
// time: current simulation time in SI units
{
    // calculate p[k][i], the probability that ion 'i' occupies state 'k'
    vector<vector<double>> p(numStates);
    for (int k = 0; k < numStates; k++){ // for each basis state
        p[k].resize(N0,0.); // initialize population to zero
    }
    for (int k = 0; k < numStates; k++){ // for each basis state
        for (int i = 0; i < N0; i++){ // for each ion
            p[k][i] = pow(norm(bwvFn[k].t()*wvFn[i]),2.); // get probability to be in state |k>
        }
    }

    // bin the state populations into pBin[k][j], the mean state population of state |k> for velocity bin 'j'
    int numBins{static_cast<int>(vBin.size())};
    vector<vector<double>> pBin(numStates);
    for (int k = 0; k < numStates; k++){
        pBin[k].resize(numBins,0.);
    }

    // compute number of particles in each bins
    double dV{vBin[1] - vBin[0]};
    bool cond1, cond2;
    for (int i = 0; i < vBin.size(); i++){ // for each velocity bin
        for (int k = 0; k < numStates; k++){ // for each basis state
            for (int j = 0; j < N0; j++){ // for each particle
                cond1 = vx[j] < (vBin[i]+dV/2.); // check if particle x-velocity is within current bin
                cond2 = vx[j] > (vBin[i]-dV/2.); // convert particle velocity to SI for comparison
                if (cond1 && cond2){
                    pBin[k][i] += p[k][j]/N0;
                }
            }
        }
    }

    // append pBin[k][j], the mean state populations of state |k> within bin 'j' to file
    string delim{","};
    ofstream outFile(filePath,ofstream::app); // append to file
    for (int k = 0; k < numStates; k++){ // for each basis state
        outFile << time << delim;
        for (int j = 0; j < numBins; j++){ // for each velocity bin
            if (j < numBins-1){
                outFile << pBin[k][j] << delim;
            }
            else{
                outFile << pBin[k][j] << endl;
            }
        }
    }

}
