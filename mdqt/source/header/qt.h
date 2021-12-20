#ifndef QT_H
#define QT_H

#include "constants.h"          // import cts:: namespace
#include "wignerSymbols.h"      // imports 'WignerSymbols::wigner3j'
#include "utils.hpp"            // imports user-defined utility functions
#include <math.h>               // imports 'pow' and 'M_PI'
#include <vector>               // imports std::vector
#include <string>               // imports std::string capability
#include <random>               // imports std::default_random_engine
#include <iostream>             // imports std::cout and std::endl
#include <armadillo>            // imports arma::cx_mat

// declare namespaces
using namespace std;
using namespace WignerSymbols;
using namespace arma;

// begin definition of QT class
class QT
{
    public:
        const int numStates{5}; // total number of basis states -> must match length of n, l, s, j, m
        const double mI{cts::mI}; // Sr+ mass in kg

        //*** NATURAL UNITS FOR QT CODE ***//
        const double timeUnit{1./cts::gam422};            // units: s
        const double freqUnit{cts::gam422};               // units: s^{-1}
        const double energyUnit{cts::hbar*cts::gam422};   // units: J
        const double lengthUnit{cts::lam422/(2.*cts::pi)}; // units: m
        const double velUnit{lengthUnit/timeUnit};        // units: m/s

        //*** PARAMETERS INITIALIZED WITH CLASS ***//
        const int N0;       // number of ions in simulation
        const double Om;    // Rabi frequency of imaging laser with units freqUnit
        const double theta; // angle that projection of B in x-y plane subtends from y-axis of lab frame (radians)
        const double phi;   // angle that B subtends from x-y plane (radians)
        const double B;     // magnitude of magnetic field (tesla)
        const double det;   // detuning of 422 imaging laser from unperturbed resonance with units freqUnit
        const string pol;   // LIF-laser polarization option: <lin>, <left>, or <right>
        const double P0;    // electron-spin polarization of LIF ground state 

        vector<cx_mat> wvFn; // wvFn[i] is wave function for particle i (column vector)

        //*** QUANTUM NUMBERS FOR EACH BASIS STATE ***//
        const vector<double> n{5., 5., 5., 5., 4.}; // principal quantum number
        const vector<double> l{0., 0., 1., 1., 2.}; // orbital angular momentum
        const vector<double> s{0.5, 0.5, 0.5, 0.5, 0.5}; // spin angular momentum
        const vector<double> j{0.5, 0.5, 0.5, 0.5, 1.5}; // total angular momentum: J=L+S
        const vector<double> m{-0.5, 0.5, -0.5, 0.5, 1.5}; // projection quantum number for J
        const vector<int> manifold{1, 1, 2, 2, 3}; // unique integer to label which manifold each state belongs to

        //*** NON-CONST CLASS MEMBERS ***//
        vector<cx_mat> bwvFn; // bwvFn[i] contains a column vector for the basis state |i>
        mat J; // J[i][k] = 3-j symbol for states |i> and |k>, note J[i][k]=J[k][i]
        cx_mat eps; // eps(i,j) is laser polarization that couples state |j> to state |i>

        vector<vector<int>> decayInd; // decayInd[k] indexes decay from excited state decayInd[k][1] to ground state decayInd[k][0]
        vector<double> lam; // lam[k] is the wavelength associated with decayInd[k]
        vector<string> label; // label[k] is a unique label for decayInd[k]

        vector<double> gamMan; // gamMan[k] is decay rate from upper to lower manifold for decayInd[k] with units freqUnit
        vector<double> gam; // gam[k] is decay rate with units freqUnit corresponding to decayInd[k] for specific transition

        cx_mat decayMat; // decayMat: decay matrix of form sum_{k} \gamma_{k}c^{\dagger}_{k}*c_{k}

  // CLASS CONSTRUCTOR
    QT(double B) : QT(1,1,1,1,B,1,"1",1) {};

    QT(int N0,double Om,double theta,double phi,double B,double det,string pol,double P0):
    N0{N0}, Om{Om}, theta{theta}, phi{phi}, B{B}, det{det}, pol{pol}, P0{P0}
    {
        // define basis wavefunctions
            mat ident{mat(numStates,numStates,fill::eye)};
            bwvFn.resize(numStates);
            for (int i = 0; i < numStates; i++){
                bwvFn[i] = cx_mat(ident.col(i),mat(numStates,1,fill::zeros));
            }

        // get 3-J symbols
        J = mat(numStates,numStates,fill::zeros);
        for (int i = 0; i < numStates; i++){  // calculate elements of J(i,k) based on quantum numbers
            for (int k = 0; k < numStates; k++){
                J(i,k) = wigner3j(j[i],1.,j[k],m[i],m[k]-m[i],-m[k]);
            }
        }

        // define decayInd[i][l] where |i> is ground state, |l> is excited state
        int numDecays = 6;
        decayInd.resize(numDecays);
        gamMan.resize(numDecays,0.);
        gam.resize(numDecays,0.);
        lam.resize(numDecays,0.);
        label.resize(4);
        decayInd[0] = {0,2}; gamMan[0] = 1.; lam[0] = cts::lam422; label[0] = "pi-";
        decayInd[1] = {0,3}; gamMan[1] = 1.; lam[1] = cts::lam422; label[1] = "sig+";// sigma plus
        decayInd[2] = {1,2}; gamMan[2] = 1.; lam[2] = cts::lam422; label[2] = "sig-";// sigma minus
        decayInd[3] = {1,3}; gamMan[3] = 1.; lam[3] = cts::lam422; label[3] = "pi+";// pi plus
        decayInd[4] = {4,2}; gamMan[4] = cts::gam1092/cts::gam422; lam[4] = cts::lam1092; // to D-state
        decayInd[5] = {4,3}; gamMan[5] = cts::gam1092/cts::gam422; lam[5] = cts::lam1092; // to D-state

        // define decay matrix: \gamma_{k}*c^{\dagger}_{k}*c_{k} from paper
        mat zeromat{mat(numStates,numStates,fill::zeros)};
        decayMat = cx_mat(zeromat,zeromat);
        int gndInd;
        int excInd;
        for (int i1 = 0; i1 < decayInd.size(); i1++){
            gndInd = decayInd[i1][0];
            excInd = decayInd[i1][1];
            if (manifold[gndInd] == 1){
                gam[i1] = (2.*j[excInd]+1)*pow(J(gndInd,excInd),2.)*gamMan[i1];
            }
            else if (manifold[gndInd] == 3){
                gam[i1] = gamMan[i1];
            }
            decayMat += gam[i1]*bwvFn[excInd]*bwvFn[excInd].t();
        }

        // define laser polarization
        complex<double> epsM{}, eps0{}, epsP{};
        getLaserPolarization(epsM,eps0,epsP);
        eps = cx_mat(zeromat,zeromat);
        int dm;
        for (int i = 0; i < numStates; i++){
            for (int k = 0; k < numStates; k++){
                if ((manifold[i] != 3) && (manifold[k] != 3) && (i != k)){
                    dm = m[i] - m[k];
                    if (dm == -1){
                        eps(i,k) = epsM;
                    }
                    else if (dm == 0){
                        eps(i,k) = eps0;
                    }
                    else if (dm == 1){
                        eps(i,k) = epsP;
                    }
                }
            }
        }
    }

    //*** NON-CONST MEMBER FUNCTIONS ***//
        void initWaveFunctions(default_random_engine &generator);
        void qtstep(const double dt,const vector<double> &vx,default_random_engine &generator,vector<double> &dv);
        void evolveWaveFunction(const double dt,const int ion,const cx_mat &Heff);
        
    //*** CONST MEMBER FUNCTIONS
        void getLaserPolarization(complex<double> &epsM,complex<double> &eps0,complex<double> &epsP) const;
        void getEffectiveHamiltonian(double vx,cx_mat &Heff) const;
        double getVelocityKickFromQTForce(double dt,cx_mat &p) const;
        complex<double> getRabiCoupling(int excInd,int gndInd) const;
        double getResonantFreq(const std::string trans) const;

        void writeStatePopulations(std::filesystem::path savePath,const vector<double> &vx) const;
        void writeQTParams(std::filesystem::path savePath,double dtQT) const;
        void writeBasisStateInfo(std::filesystem::path savePath) const;
        void writeEnsStatePops(std::filesystem::path savePath,double time) const;
        void writeBinnedStatePopsVsVel(std::filesystem::path filePath,vector<double> vBin, vector<double> vx,double time) const;
        void writeEnsAvgDensityMatrix(std::filesystem::path savePath,double time) const;
        
    };

#endif // QT_H
