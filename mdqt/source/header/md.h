#ifndef MD_H
#define MD_H

// The "Md" implements the molecular-dynamics part of the MDQT algorithm

#include "utils.hpp"    // user-defined utility functions
#include "constants.h"  // import constants contained in cts::
#include <vector>       // imports std::vector
#include <string>       // imports std::string capability
#include <math.h>       // imports 'pow' and 'M_PI'
#include <random>       // imports std::default_random_engine
#include <iostream>     // std::cout
#include <armadillo>    // for complex matrix math
#include <filesystem>   // std::filesystem

using namespace std;

class Md
{
    public: // all class members and member functions will be public
        // class members required as inputs
            const int N0;       // number of particles in simulation box
            const double n;     // plasma number density (units: m^-3)
            const double Ge;    // electron Coulomb coupling parameter

        // class members determined by class inputs
            const double L;     // length of simulation box (units: a) as defined by minimum image convention
            const double wPE;   // einstein frequency (units: s^-1)
            const double a;     // plasma interparticle spacing (units: m) (e.g. Wigner-Seitz radius)
            const double Ec;    // natural unit of energy (units: J), this is average nearest neighbor Coulomb energy
            const double Te;    // electron temperature in K
            const double lDeb;  // debye screening length (units: a)
            const double kap;   // plasma screening parameter kappa (dimensionless)

        // class members not initialized with class
            vector<vector<double>> R;     // R[i][j] holds ion positions (a) for dimension 'i' and ion 'j'
            vector<vector<double>> V;     // V[i][j] holds ion velocities (a*wPE) for dimension 'i' and ion 'j'

        // Ion class constructor
            Md(int i1,double i2,double i3):
              N0{i1},
              L{getBoxSizeForMIC(N0)},
              n{i2},
              wPE{sqrt(n*cts::e*cts::e/(3.*cts::mI*cts::eps0))},
              a{pow(3./(4.*cts::pi*n),1./3.)},
              Ec{cts::e*cts::e/(4.*cts::pi*cts::eps0*a)},
              Ge{i3},
              Te{cts::e*cts::e/(4.*cts::pi*cts::eps0*a*cts::kB*Ge)},
              lDeb{sqrt(cts::eps0*cts::kB*Te/(n*cts::e*cts::e))/a},
              kap{1./lDeb}
            {
                //initialize vectors to zero with proper size
                R.resize(3); V.resize(3); // each has three spatial dimensions
                for (int i = 0; i < R.size(); i++){ // for each spatial dimension
                    R[i].resize(N0); V[i].resize(N0);
                    std::fill(R[i].begin(),R[i].end(),0.);
                    std::fill(V[i].begin(),V[i].end(),0.);
                }

            }

        // define member functions that modify class members
            double getBoxSizeForMIC(int N0);
            void initIonPositions(default_random_engine &generator,filesystem::path path,string option = "default");
            void initIonVelocities(default_random_engine &generator,filesystem::path path,string option = "default", double Ti = 0.5);
            void andersenThermostat(default_random_engine &generator,double dt,double nu, double Ti);
            void leapfrog(double timestep, vector<vector<double>> &dV,bool turnOffMD = false);
            void applyPeriodicBC(void);

        // define functions that extract information from intrinsic plasma properties
            void getYukawaForce(vector<vector<double>> &F) const;
            double getTotalYukawaPotentialEnergy(void) const;
            vector<double> getKineticEnergy(void) const;

        // define save functions
            vector<double> getAndWriteVelocityBins(string savePath,double Ti,double Tdih) const;
            void writePlasmaParameters(std::filesystem::path filePath,double dtMD) const;
            void writeIonPosAndVel(std::filesystem::path savePath) const;
            void writePlasmaEnergy(std::filesystem::path savePath, std::filesystem::path fileName, double time, vector<double> Ekin, double Epot) const;
            void writeVelocityDistribution(std::filesystem::path savePath,vector<double> vBin, double time) const;
};

#endif
