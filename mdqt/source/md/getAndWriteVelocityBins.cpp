#include <md.h>
using namespace std;
// bin ion velocity distribution and write to file
vector<double> Md::getAndWriteVelocityBins(string savePath,double Ti,double Tdih) const
// Ti: user-specified ion temperature
// Tdih: expected disorder-induced heating temperature based on initial density
// return: velocity bins with SI units
{
    // use highest expected ion temperature to assign velocity bins for entire simulation.
    // if the specified plasma temperature is higher than the DIH temp, use that.
    // if the DIH temperature is higher than the initial ion temp, use that instead.
    double T{};
    if (Ti > Tdih) T = Ti;
    else T = Tdih;

    double sigV{sqrt(cts::kB*T/cts::mI)}; // thermal width of velocity distribution with SI units

    // using sigV, define vBin[j], the velocity in m/s of bin 'j'.
    // note that I use an odd number of bins so that I include v = 0 as a bin
    double binMax{3*sigV};
    double binMin{-3*sigV};
    int numBins{51};
    vector<double> vBin{linspace<double>(binMin,binMax,numBins)};

    // write vBin[j], the velocity bins in m/s to csv file
    string delim{","};
    ofstream outFile(savePath);
    for (int j = 0; j < numBins; j++){ // for each velocity bin, write velocity to file
        if (j < numBins-1) outFile << vBin[j] << delim;
        else outFile << vBin[j] << endl;
    }
    return vBin;
}
