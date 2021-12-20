#include "qt.h"

using namespace std;

// this function saves the basis-state-related information (i.e., the quantum numbers)
// rows correspond to basis state
// columns are organized as [n l s j m manifold]

void QT::writeBasisStateInfo(std::filesystem::path savePath) const
{
    // define file delimiter
    string delim{","};

    // open file stream, this overwrites pre-existing files
    string fileName{"basisStates.csv"};
    ofstream outFile(savePath/fileName);

    // write information for each of 'i' basis states for each quantity q[i],
    // which are vectors indexed by 'i'
    for (int i = 0; i < numStates; i++){
        outFile << n[i] << delim;
        outFile << l[i] << delim;
        outFile << s[i] << delim;
        outFile << j[i] << delim;
        outFile << m[i] << delim;
        outFile << manifold[i] << endl;
    }
    outFile.close();
}
