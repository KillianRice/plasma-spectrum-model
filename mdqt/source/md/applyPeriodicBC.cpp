#include <md.h>
// ensure all particle positions lie within simulation box
void Md::applyPeriodicBC(void)
{
    for (int i = 0; i < R[0].size(); i++){ // for each particle, return particle to simulation box if it left
        if (R[0][i] < 0) { R[0][i] += L; }
        if (R[0][i] > L) { R[0][i] -= L; }
        if (R[1][i] < 0) { R[1][i] += L; }
        if (R[1][i] > L) { R[1][i] -= L; }
        if (R[2][i] < 0) { R[2][i] += L; }
        if (R[2][i] > L) { R[2][i] -= L; }
    }
}
