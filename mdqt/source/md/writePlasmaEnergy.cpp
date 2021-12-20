#include <md.h>

using namespace std;

// This function saves the kinetic and potential energy per particle of the plasma. Each row corresponds to a given time,
// and the columns are structured as: [KEx KEy KEz Epot]

void Md::writePlasmaEnergy(std::filesystem::path savePath, std::filesystem::path fileName, double time, vector<double> Ekin, double Epot) const
{
    string delim{","};
    ofstream outFile(savePath/=fileName,ofstream::app);
    outFile << time << delim;
    outFile << Ekin[0] << delim;
    outFile << Ekin[1] << delim;
    outFile << Ekin[2] << delim;
    outFile << Epot << delim;
    outFile << Ekin[1]+Ekin[0]+Ekin[2]+Epot << endl;
}
