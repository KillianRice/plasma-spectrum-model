#include <md.h>
using namespace std;
// return Yukawa force on each particle due to all other particles
void Md::getYukawaForce(vector<vector<double>> &F) const
{
    // calculate ion forces
    F.resize(3);
    for (int i = 0; i < F.size(); i++){
        F[i].resize(N0,0.);
    }
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N0-1; i++){ // for each particle
            for (int j = i+1; j < N0; j++){ // for each other particle
                double Rd[3]{}; // intiate position difference vector (units: a)
                double RdMag{}; // magnitude of relative positions
                double f{}; // temporary force variable

                for (int k = 0; k < 3; k++){ // for each spatial dimension
                    Rd[k] = R[k][i] - R[k][j]; // calculate difference in position on that axis
                    Rd[k] -= L*round(Rd[k]/L); // apply MIC to find nearest charge
                    RdMag += pow(Rd[k],2.); // sum the square of each Rd component
                }
                RdMag = pow(RdMag,1./2.); // calculate magnitude of Rd vector

                if (RdMag < L/2. && RdMag > 0.){ // calculate force if MIC charge was found
                    f = (1./RdMag+1./lDeb)*exp(-RdMag/lDeb)/pow(RdMag,2.);
                    for (int k = 0; k < 3; k++){ // for each spatial dimension, apply equal and opposite forces
                        #pragma omp atomic update  
                        F[k][i] += Rd[k]*f;
                        
                        #pragma omp atomic update
                        F[k][j] -= Rd[k]*f;
                    }
                }
            }
        }
    }
}
