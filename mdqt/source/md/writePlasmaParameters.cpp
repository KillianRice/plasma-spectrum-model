#include <md.h>
using namespace std;
// write MD parameters to file
void Md::writePlasmaParameters(std::filesystem::path savePath,double dtMD) const
// savePath: directory path (relative to 'runFil') that specifies where files should be saved
// dtMD: MD time step with units wPE
{
    string delim{","};
    std::string fileName{"plasmaParams.csv"};
    std::ofstream outFile(savePath/=fileName);
    outFile << "N0" << delim << N0 << std::endl;
    outFile << "n" << delim << n << std::endl;
    outFile << "Ge:" << delim << Ge << std::endl;
    outFile << "L:" << delim << L << std::endl;
    outFile << "wPE:" << delim << wPE << std::endl;
    outFile << "a:" << delim << a << std::endl;
    outFile << "Ec:" << delim << Ec << std::endl;
    outFile << "Te:" << delim << Te << std::endl;
    outFile << "lDeb:" << delim << lDeb << std::endl;
    outFile << "kap:" << delim << kap << std::endl;
    outFile << "dtMD" << delim << dtMD << endl;
    outFile.close();
}
