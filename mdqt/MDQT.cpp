#include "md.h"         // imports Ion class - holds classical properties of the ions
#include "qt.h"         // imports QT class - holes quantum properties of the ions
#include "constants.h"  // imports 'cts' namespace with useful constants
#include "PlasmaSettings.hpp"   // imports functionality for .settings file
#include "utils.hpp"    // imports various user-defined functions
#include <algorithm>    // imports std::fill for vectors and matrices
#include <armadillo>    // imports 'arma' namespace, allows for vector/matrix multiplication
#include <filesystem>   // imports std::filesystem

using namespace std;
using namespace arma;
using namespace WignerSymbols;

// main function - runs MDQT simulation for unique parameter set
void MDQT(PlasmaSettings pms,int array,std::filesystem::path savePath)
{
    //*** INITIALIZE RANDOM NUMBER GENERATOR ***//
    
    default_random_engine generator;   // set up default random number generator
    int currTime{static_cast<int>(time(NULL))};         // get current time as integer
    generator.seed(array+currTime);      // seed generator with current time and array number

    //*** EXTRACT INPUTS FROM .SETTINGS FILE ***//

    const double Om{pms.getvar("Om")}; // Rabi frequency of 422 nm imaging laser with units \gamma
    const double theta{pms.getvar("theta")*(cts::pi/180.)}; // angle (radians) between projection of local magnetic field in x-y plane and y axis
    const double phi{pms.getvar("phi")*cts::pi/180.}; // angle (radians) between local magnetic field and x-y plane
    const double B{pms.getvar("B")}; // magnetic field (tesla) magnitude
    const double n{pms.getvar("n")}; // plasma number density (units: m^-3)
    const double Ti{pms.getvar("Ti")}; // initial plasma temperature (units: K)
    const int N0{static_cast<int>(pms.getvar("N"))}; // particle number, note N^{1/3} = integer when initializing positions in lattice
    const double tmax{pms.getvar("tmax")}; // time (s) at which simulation ends
    const double P{pms.getvar("P")}; // electron-spin polarization of ions in ground state

    double det; // LIF-laser detuning from unperturbed resonance
    const string det_opt{pms.getopt("det_opt")}; // option for setting laser detuning
    if (strcmp(det_opt,"input")) det = pms.getvar("det_input"); // use input value from .settings
    else if (strcmp(det_opt,"trans")){ // choose det to be resonant with one of the LIF transitions
        QT temp(B);
        det = temp.getResonantFreq(pms.getopt("det_trans"));
    }

    double Ge; // electron Coulomb coupling parameter
    const string geOpt{pms.getopt("GeOrTe")}; // option for specifying Ge
    if (strcmp(geOpt,"Ge")) Ge = pms.getvar("Ge"); // either use .settings value for Ge
    else if (strcmp(geOpt,"Te")) Ge = coulomb_coupling(n,pms.getvar("Te")); // or derive Ge from inputted Te value

    const string vOption{pms.getopt("velOpt")}; // initial velocity distribution is <maxwellian> with temp Ti or <load> from another data set
    const string pOption{pms.getopt("posOpt")}; // initialize particles with <random> or <lattice> positions or <load> from another set
    const bool thermostatIons{str2bool(pms.getopt("thermostat"))}; // <true> use Andersen thermostat to maintain ion temp at Ti <false> no bath
    const double nu{pms.getvar("nu_therm")};  // collision frequency with Andersen bath (units: \omega^{pE})
    const bool turnOffMD{str2bool(pms.getopt("turnOffMD"))};   // <true> set MD force to zero <false> calculate forces as normal
    const bool recordForMD{str2bool(pms.getopt("recordForMD"))};    // <false> record depending on QT step <true> record based on sampleFreq
    const bool turnOffQT{str2bool(pms.getopt("turnOffQT"))};   // <true> do not use QT algorithm to propagate wave functions
    const string pol{pms.getopt("pol")}; // LIF-laser polarization <lin>, <left>, or <right>

    // define time steps - later on I will ensure that their ratio is an integer
    double dtQT{pms.getvar("dt_QT")};   // QT time step with units \gamma^{-1}
    double dtMD{pms.getvar("dt_MD")};   // MD time step with units \omega_PE^{-1}

    // set time step for when to record data - based on either MD or QT timescale
    double temp = 200*dtMD/0.001*pow(n/2e14,1./2.);
    if (temp < 1.) temp = 1.1;
    int sampleFreq{static_cast<int>(temp)};              // number of time steps between recording data
    int sampleFreqQT{100};   // number of qt time steps in between recording data

    //*** INITIALIZE SAVE FILES ***//
        std::filesystem::path energyFileName{"energies.csv"};
        std::filesystem::path velBinsFileName{"velBins.csv"};
        std::filesystem::path velDistFileName{"velDist.csv"};
        std::filesystem::path ensStatePop{"ensAvgStatePops.csv"};
        std::filesystem::path binVelName{"binStatePopsVsVel.csv"};
        std::filesystem::path ensAvgDenMatName{"ensAvgDensityMatrix.csv"};
        {
            ofstream energyFile(savePath/energyFileName);
            energyFile.close();
            ofstream velBinsFile(savePath/velBinsFileName);
            velBinsFile.close();
            ofstream velDistFile(savePath/velDistFileName);
            velDistFile.close();
            ofstream ensPopFile(savePath/ensStatePop);
            ensPopFile.close();
            ofstream binVelFile(savePath/binVelName);
            binVelFile.close();
            ofstream ensAvgDenMatFile(savePath/ensAvgDenMatName);
            binVelFile.close();
        }

    //*** WRITE PROGRAM OPTIONS TO FILE ***//
    {
        string optsFileName{"options.csv"};
        ofstream optsFile(savePath/optsFileName);
        optsFile << "vOption" << "," << vOption << endl;
        optsFile << "pOption" << "," << pOption << endl;
        optsFile << "turnOffMD" << "," << turnOffMD << endl;
        optsFile << "turnOffQT" << "," << turnOffQT << endl;
        optsFile << "thermostatIons" << "," << thermostatIons << endl;
        optsFile << "nu_therm" << "," << nu << endl;
        optsFile << "Ti" << "," << Ti << endl;
        optsFile.close();
    }

    //*** MDQT SIMULATION ***//
        // Initialize MD class and QT class
            Md md(N0,n,Ge);
            QT qt{N0,Om,theta,phi,B,det,pol,P};

            double Tdih{2.*md.Ec/3./cts::kB*(1+md.kap/2.)};
            vector<double> vBin{md.getAndWriteVelocityBins(savePath/velBinsFileName,Ti,Tdih)};
            md.writePlasmaParameters(savePath,dtMD);
            qt.writeQTParams(savePath,dtQT);
            qt.writeBasisStateInfo(savePath);

            qt.initWaveFunctions(generator);
            std::filesystem::path loadpath{"equilCond"};
            if (strcmp(pOption,"load") || strcmp(vOption,"load")) loadpath /= getLoadPath(loadpath,Ti,Ge,n);
            md.initIonPositions(generator,loadpath,pOption);
            md.initIonVelocities(generator,loadpath,vOption,Ti);

        // Our code requires that dtMD>=dtQT and dtMD/dtQT = integer. This section alters dtQT and
        // dtMD appropriately. The new, altered values are the ones that I save, not the user inputs.
            double dtMD_SI{dtMD/md.wPE};
            double dtQT_SI{dtQT*qt.timeUnit};
            int timeStepRatio;
            if (dtMD_SI < dtQT_SI){
                dtQT = dtMD_SI/qt.timeUnit;
                timeStepRatio = 1;
            }
            else{
                timeStepRatio = static_cast<int>(dtMD_SI/dtQT_SI);
                dtMD = dtQT_SI*timeStepRatio*md.wPE;    // MD time step with units \omega_pE^{-1}
            }

        // Initialize containers for tracking changes in particle velocity
            vector<vector<double>> dvMD(3); // dvMD[i][j] holds change in velocity for dimension 'i' and particle 'j'
            vector<vector<double>> dvQT(3); // dvQT[i][j] holds change in velocity for dimension 'i' and particle 'j'
            for (int i = 0; i < dvMD.size(); i++){ // for each dimension
                dvMD[i].resize(N0,0.); // initialize as zeros
                dvQT[i].resize(N0,0.); // initialize as zeros
            }

        // Initialize other variables
            vector<double> Ekin{};  // total kinetic energy with units Ec
            double Epot{}; // total potential energy with units Ec
            double t{0.};   // time (s)
            int c0{0};  // counts the number of MD steps that have elapsed in total
            int c1{0};  // counts the number of QT steps that have elapsed in total
            int qtStepCount{timeStepRatio}; // tracks number of QT steps since last MD calculation
            vector<double> vxForQT; // ion velocity along LIF-laser propagation axis (x)

            while (t < tmax){
                // get x-velocity in SI units for use with QT code
                vxForQT.resize(N0,0.);
                for (int i = 0; i < N0; i++){
                    vxForQT[i] = md.V[0][i]*md.a*md.wPE;
                }

                // save relevant information 
                if (recordForMD){ // if basing sample frequency on MD
                    if (((c0 % sampleFreq) == 0) && (qtStepCount == timeStepRatio)){
                        Epot = md.getTotalYukawaPotentialEnergy();
                        Ekin = md.getKineticEnergy();
                        md.writePlasmaEnergy(savePath, energyFileName, t, Ekin, Epot);
                        md.writeVelocityDistribution(savePath/velDistFileName,vBin,t);
                    }
                }
                else{ // if basing sample frequency on QT
                    if ((c1 % sampleFreqQT) == 0){
                        Epot = md.getTotalYukawaPotentialEnergy();
                        Ekin = md.getKineticEnergy();
                        md.writePlasmaEnergy(savePath, energyFileName, t, Ekin, Epot);
                        qt.writeEnsStatePops(savePath/ensStatePop,t);
                        qt.writeEnsAvgDensityMatrix(savePath/ensAvgDenMatName,t);
                        qt.writeBinnedStatePopsVsVel(savePath/binVelName,vBin,vxForQT,t);
                    }
                }

                // evolve ion positions and calculate change in velocity due to position-verlet (leapfrog) algorithm
                // apply Andersen thermostat if applicable
                // the leap frog algorithm changes particle positions but does not change velocity. Instead, velocity change is recorded 
                // and applied gradually below
                    if (qtStepCount == timeStepRatio){
                        md.leapfrog(dtMD, dvMD,turnOffMD);
                        if (thermostatIons) md.andersenThermostat(generator,dtMD,nu,Ti);
                        qtStepCount = 0;
                        c0++;
                    }

                // apply QT algorithm to evolve ion wave functions and calculate the change in velocity
                // due to either the QT force or recoil from spontaneous emission
                if (!turnOffQT){
                    qt.qtstep(dtQT,vxForQT,generator,dvQT[0]);
                }

                // apply changes in velocity calculated previously
                    for (int i = 0; i < md.V.size(); i++){ // for each dimension
                        for (int j = 0; j < md.V[0].size(); j++){ // for each particle
                            md.V[i][j] += dvMD[i][j]/timeStepRatio+dvQT[i][j]/md.a/md.wPE;
                        }
                    }

                // ensure that ions are within simulation box
                    md.applyPeriodicBC();

                // recalculate important quantities and step time
                    t += dtQT*qt.timeUnit;
                    qtStepCount++;
                    c1++;
            }
    md.writeIonPosAndVel(savePath);
}
