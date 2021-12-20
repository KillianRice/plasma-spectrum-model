#include "qt.h"

// The effective Hamiltonian is stated in units of energyUnit = \hbar \gamma

void QT::getEffectiveHamiltonian(double v,cx_mat &Heff) const
// v: ion velocity with units velUnit
// Heff: effective hamiltonian with units energyUnit to be returned
{
    // define complex number i
    complex<double> i(0.,1.);

    // initialize real matrix with all zeros
    mat zeromat{mat(numStates,numStates,fill::zeros)};

    // define Zeeman Hamiltonian
    cx_mat Hz = cx_mat(zeromat,zeromat);
    double g;
    for (int i1 = 0; i1 < numStates; i1++){
        g = getLandeGFac(j[i1],l[i1],s[i1]);
        Hz(i1,i1) = cts::uB*B/energyUnit*g*m[i1];
    }

    // define atom hamiltonian including RWA shifts
    cx_mat Ha = cx_mat(zeromat,zeromat);
    Ha(2,2) = v-det;
    Ha(3,3) = v-det;
    Ha(4,4) = 2.*cts::pi*cts::c/cts::gam422*(1./cts::lam422 - 1./cts::lam1092);

    // define decay Hamiltonian, H_dec
    cx_mat H_dec = -i/2.*decayMat;

    // define dipole hamiltonian
    cx_mat Hd = cx_mat(zeromat,zeromat);
    for (int i1 = 0; i1 < numStates; i1++){
        for (int i2 = 0; i2 < numStates; i2++){
            if ((manifold[i1] != 3) && (manifold[i2] != 3)){
                Hd(i1,i2) = -getRabiCoupling(i1,i2)/2.;
            }
        }
    }

    // add real and imaginary parts
    Heff = Hz+Ha+Hd+H_dec;

}
