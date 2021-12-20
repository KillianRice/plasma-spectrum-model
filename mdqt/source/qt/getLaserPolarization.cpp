#include "qt.h"

// This function outputs the spherical tensor components of the LIF-laser polarization
// when projected into the local quantization axis set by the quadrupole magnetic field model.

using namespace std;

void QT::getLaserPolarization(complex<double> &epsM,complex<double> &eps0,complex<double> &epsP) const
// Note: 'option' allows the user to select the laser's polarization as defined in the lab frame, and is not
//          directory indicative of the laser polarization as expressed in local magnetic field basis.

// eta: dimensionless parameter ranging from 0 to 1, sets ratio between y- and z-components of laser polarization.
//      (eta = 1 --> y-polarization only) and (eta = 0 --> z-polarization only)
//      (eta = 1/sqrt(2) -> circular polarization)
// xi: phase shift of y-component of polarization in radians.
//      (xi = -pi/2 -> left-handed polarization)
//      (xi = pi/2 -> right-handed polarization)
// This function also uses the 'theta' parameter from the QT class, which specifies the angle (in radians) that
// the y-axis of the lab-frame coordinate system makes with the local magnetic field vector.

// eps[i]: i = 0
{
    // choose eta and xi based on user-defined 'option'
    double eta;
    double xi;
    if (pol.compare("lin") == 0 && B != 0){
        eta = 1.;
        xi = 0.;
    }
    else if (pol.compare("right") == 0 && B != 0){
        eta = 1./sqrt(2.);
        xi = -cts::pi/2.;
    }
    else if (pol.compare("left") == 0 && B != 0){
        eta = 1./sqrt(2.);
        xi = cts::pi/2.;
    }
    
    // define laser polarization in spherical tensor basis
    // eps = [eps0 epsP epsM]
    complex<double> i(0.,1.);
    eps0 = eta*cos(theta)*cos(phi)+exp(i*xi)*sqrt(1.-pow(eta,2.))*sin(phi);
    epsP = 1./sqrt(2.)*(-exp(i*xi)*sqrt(1.-pow(eta,2.))*cos(phi)-i*eta*sin(theta)+eta*cos(theta)*sin(phi));
    epsM = 1./sqrt(2.)*(+exp(i*xi)*sqrt(1.-pow(eta,2.))*cos(phi)-i*eta*sin(theta)-eta*cos(theta)*sin(phi));
    
    // if the magnetic field is zero, then the quantization axis is set by the LIF-laser
    if (B==0){
        if (pol.compare("lin")==0){
            epsP = 0;
            epsM = 0;
            eps0 = 1;
        }
        else if (pol.compare("left")==0){
            eps0 = 0;
            epsM = 1;
            epsP = 0;
        }
        else if (pol.compare("right")==0){
            eps0 = 0;
            epsM = 0;
            epsP = 1;
        }
    }
}
