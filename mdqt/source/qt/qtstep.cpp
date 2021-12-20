#include "qt.h"         // imports QT class
using namespace std;
using namespace arma;
// run QT step for each ion
void QT::qtstep(const double dt,const vector<double> &vx,default_random_engine &generator,vector<double> &dvx)
// dt: QT timestep with units freqUnit^-1
// vx[i]: ion velocity in SI units
// generator: random number generator seeded with current time and 'job' number,
//            to be used with distribution functions
// dvx[i]: velocity kick (m/s) due to spontaneous emission for ion 'i'
{
    uniform_real_distribution<double> distribution(0.,1.);
    vector<cx_mat> p(N0);
    vector<cx_mat> Heff(N0);
    #pragma omp parallel
    {
      #pragma omp for
      for (int i = 0; i < N0; i++){ // for each ion

        //*** OBTAIN RANDOM NUMBERS FROM RNG ***//
        // these random numbers are between [0,1] and are used for various things
        double r[3];
        #pragma omp critical
        {
          for (int j = 0; j < 3; j++){
            r[j] = distribution(generator);
          }
        }


        //*** DETERMINE WHETHER THE CURRENT ION DECAYS DURING THIS TIME STEP ***//
          // calculate dP[k] the probability that ion 'i' will decay via decay path 'k'
          // also calculate dPtot = sum_{k}dP[k]
          double dP[decayInd.size()];
          double dPtot{0.};
          int gndInd;
          int excInd;

          for (int k = 0; k < decayInd.size(); k++){ // for each allowed decay path
              excInd = decayInd[k][1];
              dP[k] = dt*norm(wvFn[i].t()*bwvFn[excInd]*bwvFn[excInd].t()*wvFn[i]*gam[k]);
              dPtot += dP[k];
          }

          // randomly decide whether a decay happens (e.g., set doesIonDecay = true/false)
          bool doesIonDecay{false};
          if (r[0] < dPtot){
            doesIonDecay = true;
          }

        //*** HANDLE WHEN ION DECAYS ***//
          if (doesIonDecay){
            // decide which decay path is used
              int decayPath{1000};  // initialize with unrealistic value for future comparison
              double ind{0.};
              for (int k = 0; k < decayInd.size(); k++){
                ind += dP[k]/dPtot;
                if (r[1] < ind){
                  decayPath = k;
                  break;
                }
              }

            // obtain dvx[i], the velocity kick for ion 'i', from photon recoil with random direction for decay path 'k'
            double sign{1.};
            if (r[2] < 0.5){
                sign = -1.;
            }
            dvx[i] = sign*cts::hbar*(2.*cts::pi/lam[decayPath]); // this is in SI units intentionally

            // set wave function of ion 'i' equal to ground state corresponding to decay path 'k'
            gndInd = decayInd[decayPath][0];
            wvFn[i] = bwvFn[gndInd];  // set real part of groundInd to one
          }

        //*** HANDLE WHEN ION DOES NOT DECAY ***//
          if (!doesIonDecay){
            p[i] = wvFn[i]*wvFn[i].t(); // form density matrix
            dvx[i] = getVelocityKickFromQTForce(dt,p[i])*velUnit; // get velocity kick from QT force in SI units
            getEffectiveHamiltonian(vx[i]/velUnit,Heff[i]); // get effective Hamiltonian
            evolveWaveFunction(dt,i,Heff[i]); // evolve wave functions according to effective Hamiltonian
          }
      }
    }
}
