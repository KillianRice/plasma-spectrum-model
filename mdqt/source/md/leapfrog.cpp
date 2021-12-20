#include <md.h>
using namespace std;

// This function implements a position-verlet leapfrog algorithm to detect the change in ion position
// and velocity due to the screened Yukawa force. The leapfrog algorithm is only valid for spatially-dependent forces. Use
// care when applying it to other potentials.

// implement leap-frog algorithm to change particle positions and record changes in velocity
void Md::leapfrog(double timestep, vector<vector<double>> &dV,bool turnOffMD)
{
    //*** FIRST HALF STEP IN POSITION ***//
        for (int i = 0; i < R.size(); i++){ // for each spatial dimension
            for (int j = 0; j < R[0].size(); j++){ // for each particle
                R[i][j] += V[i][j]*timestep/2; // update particle positions
            }
        }
        applyPeriodicBC(); // ensure particles stay within simulation box

    //*** FULL STEP IN VELOCITY
    
        // calculate forces between ions at intermediate positions
        vector<vector<double>> F(3);
        for (int i = 0; i < F.size(); i++){
        F[i].resize(N0,0.);
        }
        if (!turnOffMD) getYukawaForce(F);
        
        // calculate change in velocity due to full step
        for (int i = 0; i < R.size(); i++){ // for each spatial dimension
            for (int j = 0; j < R[0].size(); j++){ // for each particle
                dV[i][j] = F[i][j]*timestep; // record change in position due to half step
            }
        }

    //*** SECOND HALF STEP IN POSITION
        for (int i = 0; i < R.size(); i++){ // for each spatial dimension
            for (int j = 0; j < R[0].size(); j++){ // for each particle
                R[i][j] += (V[i][j]+dV[i][j])*timestep/2; // update particle positions
            }
        }
        applyPeriodicBC(); // ensure particles are within simulation box
}
