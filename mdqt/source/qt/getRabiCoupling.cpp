#include "qt.h"
// get Rabi coupling between state gndInd and excInd
complex<double> QT::getRabiCoupling(int excInd,int gndInd) const
{
    double j_gt;
    if (j[excInd] > j[gndInd]){
        j_gt = j[excInd];
    }
    else {
        j_gt = j[gndInd];
    }

    return Om*pow(-1,j_gt+j[excInd]+j[gndInd]-m[excInd])*J(excInd,gndInd)*sqrt(2.*j[excInd]+1.)*conj(eps(excInd,gndInd));
}