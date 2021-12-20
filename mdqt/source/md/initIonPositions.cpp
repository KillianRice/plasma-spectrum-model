#include <md.h>
using namespace std;
// initialize particle positions within simulation box
void Md::initIonPositions(default_random_engine &generator,std::filesystem::path path,string option)
{
    //*** INITIALIZE POSITIONS RANDOMLY ***//
    if (option.compare("random") == 0){
        // Initialize ion positions (units of a) with box of length L using pseudo-random number generator. In order to
        // ensure that ion positions are truly random, we sample ion positions within a much bigger box centered around
        // the simulation box. If the sampled position lies within the simulation box, we record the position.
            uniform_real_distribution<double> distribution (-0.5*L,1.5*L);

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < N0; i++){
                double Rs[3]{}; // current sampled position
                bool recorded{false}; // whether or not this particle's position has been recorded
                while (!recorded){
                    // sample positions
                    #pragma omp critical // generator cannot be accessed by different threats simultaneously
                    {
                        Rs[0] = distribution(generator);
                        Rs[1] = distribution(generator);
                        Rs[2] = distribution(generator);
                    }

                    // record position for this particle if Rs lies in the box
                    if (Rs[0]<=L && Rs[1]<=L && Rs[2]<=L && Rs[0]>0 && Rs[1]>0 && Rs[2]>0){
                        R[0][i] = Rs[0];
                        R[1][i] = Rs[1];
                        R[2][i] = Rs[2];

                        recorded = true;
                    }
                }
            }

        }
    }
    //*** INITIALIZE IONS WITHIN A LATTICE ***//
    else if (option.compare("lattice") == 0){
        // This section initializes the ion positions inside a 3D box occupying r = [0,L] on each axis.
        int numPos{round2int(pow(N0,1./3.))};
        double dx{L/(2*(numPos+1)+2)};
        vector<double> xopt{linspace<double>(dx,L-dx,numPos)};

        vector<vector<double>> pos;
        for (int i = 0; i < numPos; i++){
            for (int j = 0; j < numPos; j++){
                for (int k = 0; k < numPos; k++){
                    pos.push_back({xopt[i],xopt[j],xopt[k]});
                }
            }
        }

        for (int i = 0; i < N0; i++){
            for (int j = 0; j < 3; j++){
                R[j][i] = pos[i][j];
            }
        }

    }
    //*** LOAD ION POSITIONS FROM FILE ***//
    else if (option.compare("load") == 0){
        // matIn[particle][i] contains {rx ry rz vx vy vz} for each particle
        vector<vector<double>> matIn = readCSV2Double(path);

        if (matIn.size() != N0){
            cerr << "Error: MD class has a different number of particles than the equilibrated conditions that are being loaded.";
        }

        // extract the positions from matIn
        int axInd;
        for (int part = 0; part < N0; part++){
            axInd = 0;
            for (int vInd = 0; vInd < 3; vInd++){
                R[axInd][part] = matIn[part][vInd];
                axInd++;
            }
        }
    }
}
