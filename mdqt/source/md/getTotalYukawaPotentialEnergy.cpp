#include <md.h>
// return total potential energy between particles with units Ec - sum over all particles
double Md::getTotalYukawaPotentialEnergy(void) const
{
    double Epot{0};   // initialize potential energy per particle (units: Ec) to zero prior to loop

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N0-1; i++){ // for each particle
            for (int j = i+1; j < N0; j++){ // for each other particle
                double Rd[3]{}; // intiate position difference vector (units: a)
                double RdMag{}; // magnitude of relative positions

                for (int k = 0; k < 3; k++){ // for each spatial dimension
                    Rd[k] = R[k][i] - R[k][j]; // calculate difference in position on that axis
                    Rd[k] -= L*round(Rd[k]/L); // apply MIC to find nearest charge
                    RdMag += Rd[k]*Rd[k]; // sum the square of each Rd component
                }
                RdMag = sqrt(RdMag);               // calculate magnitude of Rd vector

                if (RdMag < L/2. && RdMag > 0.){ // calculate force if MIC charge was found
                    #pragma omp atomic update
                    Epot += exp(-RdMag/lDeb)/(RdMag)/N0;
                }
            }
        }
    }
    return Epot;
}
