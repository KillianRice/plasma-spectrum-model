#include "qt.h"

// Using dp/dt=m*dv/dt=F --> dv = F dt / m this function calculates the change in velocity
// due to the QT force in natural units.

double QT::getVelocityKickFromQTForce(double dt,cx_mat &p) const
// dt: QT time step with units 'timeUnit'
// p: density matrix of particle, access element (i,j) using 'p(i,j)'
// dv: velocity kick (velUnit) due to QT force
{
    double dv = 0;
    complex<double> Om_eg;
    complex<double> i(0.,1.);
    double k = 2.*cts::pi/cts::lam422;
    for (int exc = 0; exc < numStates; exc++){
        for (int gnd = 0; gnd < numStates; gnd++){
            if ((manifold[exc]==2) && manifold[gnd]==1 ){
                Om_eg = getRabiCoupling(exc,gnd);
                dv += real(i*cts::hbar*k/2./mI/velUnit*(Om_eg*conj(p(exc,gnd))-conj(Om_eg)*p(exc,gnd))*dt);
            }
        }
    }
    return dv;
}
