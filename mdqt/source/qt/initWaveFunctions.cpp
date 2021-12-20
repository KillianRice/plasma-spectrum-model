#include "qt.h"     // imports QT class
using namespace std;
// initialize wavefunction for each particle
void QT::initWaveFunctions(default_random_engine &generator)
{
    wvFn.resize(N0); // initialize vector of wave functions to have element for each particle
    uniform_real_distribution<double> distribution (0.,1.);  // define uniform distribution between 0 and 1
    cx_mat realWvFn, imagWvFn;
    for (int i = 0; i < N0; i++){ // for each particle
        double p0 = (1.+P0)/2.; // fraction of ions that occupy m_j=+1/2 ground state
        if (distribution(generator) < p0) realWvFn = bwvFn[1];
        else realWvFn = bwvFn[0];
        imagWvFn = cx_mat(mat(numStates,1,fill::zeros),mat(numStates,1,fill::zeros));
        wvFn[i] = realWvFn+imagWvFn;
    }
}
