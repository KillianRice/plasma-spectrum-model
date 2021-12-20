#include "qt.h"

// return Zeeman-shifted resonant frequency for requested transition with units energyUnit
double QT::getResonantFreq(const std::string trans) const
{
    double dEz;
    bool found{false};
    for (int i = 0; i < label.size(); i++){
        if (strcmp(trans,label[i])){
            int exc = decayInd[i][1]; // index for excited state
            int gnd = decayInd[i][0]; // index for ground state
            double gL_exc = getLandeGFac(j[exc],l[exc],s[exc]); // lande G factor for excited state
            double gL_gnd = getLandeGFac(j[gnd],l[gnd],s[gnd]); // lande G factor for ground state
            dEz = getZeemanShift(B,gL_exc,m[exc]) - getZeemanShift(B,gL_gnd,m[gnd]);
            dEz /= energyUnit;
            found = true;   
            break;  
        }
    }
    if (!found) std::cerr << "Requested transition does not correspond to label in qt.h" << std::endl;
    return dEz;
}