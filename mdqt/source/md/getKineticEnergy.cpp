#include <md.h>
using namespace std;
// return total kinetic energy for all particles summed over all spatial dimensions
vector<double> Md::getKineticEnergy(void) const
{
    vector<double> Ek(R.size(),0.);
    for (int i = 0; i < V.size(); i++){ // for each spatial dimension
        for (int j = 0; j < V[0].size(); j++){ // for each particle
            Ek[i] += 0.5*V[i][j]*V[i][j]/N0; // add kinetic energy
        }
    }
    return Ek;
}
