#include "qt.h"

void QT::writeEnsAvgDensityMatrix(std::filesystem::path savePath,double time) const
{
    // form density matrix from particle wave functions for each ion
    // p[i] is the density matrix for particle |i>
    vector<cx_mat> p(N0);
    for (int i = 0; i < N0; i++){
        p[i] = wvFn[i]*wvFn[i].t();
    }

    // calculate the ensemble-averaged density matrix
    // pEns is a 4x4 cx_mat
    mat zeroMat{mat(numStates,numStates,fill::zeros)};  // used to initialize pEns with zeros
    cx_mat pEns{zeroMat,zeroMat};                       // initialize pEns to zeros
    for (int i = 0; i < N0; i++){                       // average the single-particle density matrices
        pEns += p[i]/N0;
    }

    // reformat pEns into 'data', a vector for writing to file
    // each row of the file will contain density matrix elements for a single time
    // each row will contain a series of real/imaginary pairs corresponding to density matrix elements
    vector<double> data{};
    data.reserve(numStates*numStates*2);
    for (int i = 0; i < numStates; i++){
        for (int j = 0; j < numStates; j++){
            data.push_back(pEns(i,j).real());
            data.push_back(pEns(i,j).imag());
        }
    }
    data.shrink_to_fit();

    // append 'data' vector to file using comma as delimiter
    string delim{","};
    ofstream outFile(savePath,ofstream::app);
    outFile << time << delim;
    for (int i = 0; i < data.size(); i++){
        if (i < data.size()-1){
            outFile << data[i] << delim;
        }
        else{
            outFile << data[i] << endl;
        }
    }

    // close file stream
    outFile.close();
}
