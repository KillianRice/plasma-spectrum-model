#include "qt.h"

using namespace arma;
using namespace std;

void QT::evolveWaveFunction(const double dt,const int ion,const cx_mat &Heff)
{
    // define complex number i = \sqrt{-1}
    complex<double> I(0., 1.);

    // define identity matrix
    mat ident{mat(numStates,numStates,fill::eye)};

    // initialize dPtot, the total probability of the ion to decay
    double dPtot;

    // first iteration of Runge-Kutta
    dPtot = dt*norm(wvFn[ion].t()*decayMat*wvFn[ion]);
    cx_mat k1{((ident-I*Heff*dt)/sqrt(1.-dPtot)-ident)*wvFn[ion]};
    cx_mat wvFn1{wvFn[ion]+k1/2.};

    // second iteration of Runge-Kutta
    dPtot = dt*norm(wvFn1.t()*decayMat*wvFn1);
    cx_mat k2{((ident-I*Heff*dt)/sqrt(1.-dPtot)-ident)*wvFn1};
    cx_mat wvFn2{wvFn[ion]+k2/2.};

    // third iteration of Runge-Kutta
    dPtot = dt*norm(wvFn2.t()*decayMat*wvFn2);
    cx_mat k3{((ident-I*Heff*dt)/sqrt(1.-dPtot)-ident)*wvFn2};
    cx_mat wvFn3{wvFn[ion]+k3};

    // fourth iteration of Runge-Kutta
    dPtot = dt*norm(wvFn3.t()*decayMat*wvFn3);
    cx_mat k4{((ident-I*Heff*dt)/sqrt(1.-dPtot)-ident)*wvFn3};

    // evolve wvFn[ion] with Runge-Kutta calculations, then normalize wave function
    wvFn[ion] += (k1+2.*k2+2.*k3+k4)/6.;
    wvFn[ion] /= sqrt(norm(wvFn[ion].t()*wvFn[ion]));

}
