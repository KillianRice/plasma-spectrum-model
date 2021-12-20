#include <md.h>
using namespace std;
// use Andersen thermostat to collisionally relax ion temperature to bath temp
void Md::andersenThermostat(default_random_engine &generator,double dt,double nu, double Ti)
{
    // set up maxwell-boltzmann generator
    double sigV{pow(cts::kB*Ti/cts::mI,1./2.)/a/wPE};   // standard deviation of velocity distribution (units: a*wPE)
    double v0{0.};                                      // center of velocity distribution (units: a*wPE)
    normal_distribution<double> maxwell_dist(v0,sigV);     // configure normal distribution

    // set up uniform distribution generator
      uniform_real_distribution<double> distribution(0.,1.);

    double P{nu*dt}; // probability that each particle will collide with bath

    double r;
    for (int i = 0; i < N0; i++){   // iterate over all paricles
        r = distribution(generator);
        if (r < P){ // if particle collides with bath
            for (int j = 0; j < 3; j++){ // for each spatial dimension
                V[j][i] = maxwell_dist(generator); // sample velocity from maxwell-boltzmann
            }
        }
    }
}
