#include <md.h>
using namespace std;

// This function writes the ion positions and velocities to a .dat file, in which each row
// corresponds to a different particle.

// Each row is structured as: [rx ry rz vx vy vz]

// All quantities recorded in the dimensionless units defined by the ion class.

void Md::writeIonPosAndVel(std::filesystem::path savePath) const
{

    string delim{","};

    string fileName{"ionPosAndVel.csv"};
    ofstream outFile(savePath/=fileName);

    for (int i = 0; i < R[0].size(); i++){ // for each particle
        outFile << R[0][i] << delim;
        outFile << R[1][i] << delim;
        outFile << R[2][i] << delim;
        outFile << V[0][i] << delim;
        outFile << V[1][i] << delim;
        outFile << V[2][i] << endl;
    }

    outFile.close();
}
