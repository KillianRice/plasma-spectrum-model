#include <md.h>
using namespace std;
// initialize particle velocities
void Md::initIonVelocities(default_random_engine &generator,filesystem::path path,string option, double Ti)
{
    //*** INITIALIZE VELOCITIES WITH MAXWELL-BOLTZMANN WITH TEMPERATURE Ti ***//
    if (strcmp(option,"maxwellian")){
        // When MD is combined with QT, we are often interested in how the ion state populations change as a function of
        // velocity. Following DIH, the ions usually span a suitable range of velocity space. However, for very dilute plasmas
        // (or when the interactions are turned off entirely), it is sometimes necessary to artificially set the ion velocities
        // to span a relevant range of velocity space. This option allows us to do that.

        // set up maxwell-boltzmann number generator
        double sigV{pow(cts::kB*Ti/cts::mI,1./2.)/a/wPE};   // standard deviation of velocity distribution (units: a*wPE)
        double v0{0.};                                      // center of velocity distribution (units: a*wPE)
        normal_distribution<double> maxwell_dist(v0,sigV);     // configure normal distribution

        // sample particle positions
        for (int i = 0; i < N0; i++){ // for each particle
            for (int j = 0; j < 3; j++){ // for each spatial dimension
                V[j][i] = maxwell_dist(generator);
            }
        }
    }
    //*** LOAD ION VELOCITIES ***//
    else if (strcmp(option,"load")){
        // matIn[part][i] contains {rx ry rz vx vy vz} for each particle
        vector<vector<double>> matIn = readCSV2Double(path);

        int axInd;
        for (int part = 0; part < N0; part++){ // for each particle
            axInd = 0;
            for (int vInd = 3; vInd < 6; vInd++){
                V[axInd][part] = matIn[part][vInd];
                axInd++;
            }
        }
    }
}
