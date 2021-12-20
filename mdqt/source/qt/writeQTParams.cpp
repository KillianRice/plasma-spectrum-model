#include "qt.h"

using namespace std;

void QT::writeQTParams(std::filesystem::path savePath,double dtQT) const
// savePath: directory path (relative to 'runFil') that specifies where files should be saved
// dtQT: QT time step with units timeUnit (should be \gamma^{-1})
{
    string fileName{"qtParams.csv"};
    ofstream outFile(savePath/fileName);
    string delim{","};
    outFile << "numStates" << delim << numStates << endl;
    outFile << "mI" << delim << mI << endl;
    outFile << "timeUnit" << delim << timeUnit << endl;
    outFile << "freqUnit" << delim << freqUnit << endl;
    outFile << "energyUnit" << delim << energyUnit << endl;
    outFile << "lengthUnit" << delim << lengthUnit << endl;
    outFile << "velUnit" << delim << velUnit << endl;
    outFile << "N0" << delim << N0 << endl;
    outFile << "Om" << delim << Om << endl;
    outFile << "theta" << delim << theta << endl;
    outFile << "B" << delim << B << endl;
    outFile << "det" << delim << det << endl;
    outFile << "dtQT" << delim << dtQT << endl;
    outFile << "pol" << delim << pol << endl;
    outFile << "P" << delim << P0 << endl;
    outFile.close();
}
