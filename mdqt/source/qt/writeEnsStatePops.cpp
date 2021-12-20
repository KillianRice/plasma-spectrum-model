#include "qt.h"

void QT::writeEnsStatePops(std::filesystem::path savePath,double time) const
// savePath: file path relative to 'runFile' (and including file name) for where
//      the ensemble-averaged ion state populations will be saved
// time: current time within the simulation in SI units
{
    // initialize data[i], where data.size() = numStates*2+1 and holds the mean/std of
    // the ensemble-averaged ion state populations. For each basis state |k>, the population mean
    // , mean[k], and standard devation std[k] will be calculated and stored in data[i] as
    // data = [mean[1],std[1],mean[2],std[2],....mean[k],std[k]]
    vector<double> data(numStates*2,0.);

    // calculate pops(i,k), the ion state population for ion 'i' and basis state |k>
    mat pops{mat(N0,numStates,fill::zeros)};
    for (int i = 0; i < N0; i++){
        for (int k = 0; k < numStates; k++){
            pops(i,k) = pow(norm(bwvFn[k].t()*wvFn[i]),2.);
        }
    }

    // get ensemble statistics and store in 'data'
    mat popMean{mean(pops,0)};  // the zero specifies that I want to average over rows
    mat popStd{stddev(pops,0)};

    // place state population information within 'data'
    int ind = 0;
    for (int k = 0; k < numStates; k++){
        data[ind] = popMean[k];
        ind++;
        data[ind] = popStd[k];
        ind++;
    }

    // specify file delimiter
    string delim{","};

    // open file stream in append mode
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
