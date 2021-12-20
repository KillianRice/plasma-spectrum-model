#include <md.h>

using namespace std;

void Md::writeVelocityDistribution(std::filesystem::path savePath,vector<double> vBin, double time) const
// savePath: directory path to
// vBin: velocity bins in SI units
// time: current simulation time in SI units
{
    // compute number of particles in each bins
    vector<double> numInBins(vBin.size(),0);
    double dV{vBin[1] - vBin[0]};
    bool cond1, cond2;
    for (int i = 0; i < vBin.size(); i++){ // for each velocity bin
        for (int j = 0; j < N0; j++){ // for each ion
            cond1 = V[0][j]*a*wPE < (vBin[i]+dV/2.); // check if particle x-velocity is within current bin
            cond2 = V[0][j]*a*wPE > (vBin[i]-dV/2.); // convert particle velocity to SI for comparison
            if (cond1 && cond2){
                numInBins[i] += 1./N0;
            }
        }
    }

    // write number in each bin to file
    string delim{","};
    ofstream outFile(savePath,ofstream::app);
    outFile << time << delim;
    for (int j = 0; j < vBin.size(); j++){
        if (j < vBin.size()-1){
            outFile << numInBins[j] << delim;
        }
        else{
            outFile << numInBins[j] << endl;
        }
    }
}
