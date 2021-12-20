#include "utils.hpp"
#include "constants.h"
#include <iostream>

bool str2bool(const std::string& str)
{
    bool result;
    if (strcmp(str,"false")) result = false;
    else if (strcmp(str,"true")) result = true;
    else std::cerr << "Input string must be <true> or <false>" << std::endl;
    return result;
}

// parse command line arguments
std::string getCommandLineArg(int argc, char* argv[], std::string short_flag, std::string long_flag)
{
  int num_args = argc - 1;
  std::vector<std::string> arguments(argv+1, argv+argc);
  std::string result;
  for(int i=0; i<num_args; i++){
    if (strcmp(arguments[i],short_flag) || strcmp(arguments[i],long_flag)){
      result = arguments[i+1];
    }
    i++;
  }
  return result;
}

// round double to int
int round2int(const double& a) { return static_cast<int>(a < 0 ? a - 0.5 : a + 0.5); }

// bin <vec_in> into <num_bins> linearly spaced bins between min and max
std::vector<double> bin_vector(const std::vector<double>& vec_in,const std::vector<double>& bins)
{
    std::vector<double> num_in_bin(bins.size(),0.);
    double bin_spacing = bins[1] - bins[0];
    
    #pragma omp parallel for
    for (int i = 0; i < bins.size(); i++){ 
        for (int j = 0; j < vec_in.size(); j++){
            bool cond1, cond2, in_bin;
            cond1 = vec_in[j] < (bins[i] + bin_spacing/2.); 
            cond2 = vec_in[j] > (bins[i] - bin_spacing/2.); 
            in_bin = cond1 && cond2;
            if (in_bin) num_in_bin[i] += 1./vec_in.size();
        }
    }
    return num_in_bin;
}

// get unique combinations of each row of mat_in and each element of vec_2
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::vector<std::string>>& mat_in,const std::vector<std::string>& vec_2)
{
    std::vector<std::vector<std::string>> mat_out(mat_in.size()*vec_2.size());
    int row = 0;
    for (int i = 0; i < mat_in.size(); i++){
        for (int j = 0; j < vec_2.size(); j++){
            mat_out[row].reserve(mat_in[i].size()+1);
            for (int k = 0; k < mat_in[i].size(); k++) mat_out[row].push_back(mat_in[i][k]);
            mat_out[row].push_back(vec_2[j]);
            row++;
        }
    }

    return mat_out;
}

// get unique combinations of the elements in the two input vectors
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::string>& vec_in,const std::vector<std::string>& vec_2)
{
    std::vector<std::vector<std::string>> mat_out(vec_in.size()*vec_2.size());
    int row = 0;
    for (int i = 0; i < vec_in.size(); i++){
        for (int j = 0; j < vec_2.size(); j++){
            mat_out[row].push_back(vec_in[i]);
            mat_out[row].push_back(vec_2[j]);
            row++;
        }
    }

    return mat_out;
}

// compute the norm (i.e., square root of a vectors inner product with itself)
double euclidean_norm(const std::vector<double>& vec_in)
{
    double norm;
    for (auto val : vec_in) norm += pow(val,2.);
    norm = sqrt(norm);
    return norm;
}

// read file with comma delimiter, "=" functions as ",", all text after "%" is ignored
std::vector<std::vector<std::string>> readCSV(std::filesystem::path filePath)
// filePath: full path to .csv file to read from, formatted as std::filesystem::path
{
    // check that file exists and is not empty
    if (!std::filesystem::exists(filePath) || std::filesystem::is_empty(filePath)){
        std::cerr << "The specified .csv file is either empty or does not exist." << std::endl;
    }

    // initialize data with reasonable buffer size
    std::vector<std::vector<std::string>> data;
    data.reserve(100);

    // open input stream to file
    std::ifstream fileStream;            // initialize empty stream object
    while(!fileStream.is_open()){             // while the file is not open...
        fileStream.open(filePath);  // try to open file
    }
    
    // read data from file line-by-line: % comments out line, read after :, spaces are ignored
    while (fileStream.good()){ // while stream is open and has no errors
        // read in current line from .csv file
        std::string currLine; // initialize container for current line of .csv file
        std::getline(fileStream,currLine); // read in current line

        // if first character is "%", line is ignored
        std::size_t pos = currLine.find("%");
        if (pos == 0) currLine.clear();
        else if (pos != std::string::npos) currLine = currLine.substr(0,pos-1);

        // convert '#' to ','
        while ((pos = currLine.find("=")) != std::string::npos) currLine.replace(pos,1,",");

        // remove spaces from string
        currLine.erase(std::remove(currLine.begin(),currLine.end(),' '),currLine.end());
        
        std::istringstream ss(currLine); // create stream to current line

        // parse .csv delimited values one at a time
        std::vector<std::string> currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){ // for each delimited value in currLine
            std::string s;
            std::getline(ss,s,','); // place delimited value in 's'
            currVec.push_back(s); // place that element in currVec

            if(currVec.size() == currVec.capacity()){ // ensure that reserved size for currVec is large enough
                currVec.reserve(2.*currVec.capacity()); 
            }
        }
        currVec.shrink_to_fit(); // remove excess space

        // store parsed line into data and ensure the reserved storage is large enough
        
        if (!currLine.empty()){
            data.push_back(currVec); // store in data
            if(data.size() == data.capacity()){ // if actual size has reached buffer size
            data.reserve(2.*data.capacity()); // double buffer size
        }
        
    }
    }

    // close file stream and shrink data buffer to actual size
    fileStream.close();
    data.shrink_to_fit();

    return data;
}

// read CSV to vector vector double
std::vector<std::vector<double>> readCSV2Double(std::filesystem::path filePath)
// filePath: full path to .csv file to read from, formatted as std::filesystem::path
{
    // check that file exists and is not empty
    if (!std::filesystem::exists(filePath) || std::filesystem::is_empty(filePath)){
        std::cerr << "The specified .csv file is either empty or does not exist." << std::endl;
    }

    // initialize data with reasonable buffer size
    std::vector<std::vector<double>> data;
    data.reserve(100);

    // open input stream to file
    std::ifstream fileStream;            // initialize empty stream object
    while(!fileStream.is_open()){             // while the file is not open...
        fileStream.open(filePath);  // try to open file
    }
    
    // read data from file line-by-line
    while (fileStream.good()){ // while stream is open and has no errors
        // read in current line from .csv file
        std::string currLine; // initialize container for current line of .csv file
        if (!std::getline(fileStream,currLine)){ // read in current line
            break;
        }
        std::istringstream ss(currLine); // create stream to current line

        // parse .csv delimited values one at a time
        std::vector<double> currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){ // for each delimited value in currLine
            std::string s;
            std::getline(ss,s,','); // place delimited value in 's'
            currVec.push_back(atof(s.c_str())); // place that element in currVec

            if(currVec.size() == currVec.capacity()){ // ensure that reserved size for currVec is large enough
                currVec.reserve(2.*currVec.capacity()); 
            }
        }
        currVec.shrink_to_fit(); // remove excess space

        // store parsed line into data and ensure the reserved storage is large enough
        data.push_back(currVec); // store in data
        if(data.size() == data.capacity()){ // if actual size has reached buffer size
        data.reserve(2.*data.capacity()); // double buffer size
    }
    }

    // close file stream and shrink data buffer to actual size
    fileStream.close();
    data.shrink_to_fit();

    return data;
}

// get path to closest set of equilibrated initial conditions
std::filesystem::path getLoadPath(std::filesystem::path path,double Ti, double Ge, double n)
{
    std::string ind;
    std::vector<std::string> files;
    enum Var {v_Ti,v_Ge,v_n,v_num};
    std::vector<std::vector<double>> var(v_num);
    for (auto i : var) i.reserve(1000);
    files.reserve(1000);
    int ind1, ind2, ind3, ind4;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        ind = entry.path();
        if (std::filesystem::is_directory(ind)){
            files.push_back(ind.substr(path.string().size()+1));
            ind1 = files.back().find("Ti");
            ind2 = files.back().find("_n");
            ind3 = files.back().find("_Ge");
            ind4 = files.back().find("_N");
            var[v_Ti].push_back(stod(files.back().substr(ind1+2,ind2-1)));
            var[v_n].push_back(stod(files.back().substr(ind2+2,ind3-1)));
            var[v_Ge].push_back(stod(files.back().substr(ind3+3,ind4-1)));

        }

    }
    files.shrink_to_fit();
    for (auto i : var) var.shrink_to_fit();

    // find which file is closest to initial conditions
    std::vector<double> diff(files.size());
    {
        for (int i = 0; i < diff.size(); i++){
            diff[i] += abs(Ti - var[v_Ti][i])/Ti;
            diff[i] += abs(n - var[v_n][i])/n;
            diff[i] += abs(Ge - var[v_Ge][i])/Ge;
        }
    }

    std::filesystem::path loadpath;
    {
        double min = std::numeric_limits<double>::max();
        int min_loc = -1;
        for (int i = 0; i < diff.size(); i++){
            if (diff[i] < min){
                min = diff[i];
                min_loc = i;
            }
        }
        return files[min_loc] + "/ionPosAndVel.csv";
    }

}

// compare <str1> and <str2>
bool strcmp(std::string str1, std::string str2)
{
    bool is_equal;
    if (str1.compare(str2) == 0) is_equal = true;
    else is_equal = false;
    return is_equal;
}

// ratio between average Coulomb interaction energy and thermal energy (inputs in SI)
double coulomb_coupling(const double& n,const double& T)
// n: plasma density (m^-3)
// T: temperature (K) of a particular species (i.e., electron or ion)
{ return nearest_coulomb_pot(n)/(cts::kB*T); }

// average interparticle spacing of gas with density n (SI units)
double wigner_seitz_radius(const double& n) 
// n: plasma density (m^-3)
{ return pow(3./(4.*cts::pi*n),1./3.); } 

// plasma oscillation frequency (rad/s, inputs SI)
double plasma_freq(const double& n,const double& m) 
// n: plasma density (m^-3)
// m: particle mass (kg)
{ return sqrt(n*pow(cts::e,2.)/m/cts::eps0); }

// einsten frequency (rad/s, inputs SI)
double einstein_freq(const double& n,const double& m)
// n: plasma density (m^-3)
// m: particle mass (kg)
{ return plasma_freq(n,m)/sqrt(3.); }

// average Coulomb potential between nearest neighbors, SI units
double nearest_coulomb_pot(const double& n) 
// n: plasma density (m^-3)
{ return pow(cts::e,2.)/(4.*cts::pi*cts::eps0*wigner_seitz_radius(n)); } 

// debye length of a particular plasma species, SI units
double debye_length(const double& n,const double& T) 
// n: plasma density (m^-3)
// T: temperature (K) of a particular plasma species
{ return sqrt(cts::eps0*cts::kB*T/(n*pow(cts::e,2.)));}

// plasma screening parameter for electrons, dimensionless, inputs SI
double screening_parameter(const double& n, const double& Te)
// n: plasma density (m^-3)
// Te: electron temperature (K)
{ return wigner_seitz_radius(n)/debye_length(n,Te); } 

// equilibrium plasma temperature due to disorder-induced heating, inputs SI
double dih_temp(const double& n,const double& Te) 
// n: plasma density (m^-3)
// Te: electron temperature (K)
{ return 2.*nearest_coulomb_pot(n)/3./cts::kB*(1.+screening_parameter(n,Te)/2.); }

// timescale for a Gaussian UCNP expansion into vacuum, SI units
double get_tau_exp(const double& sig,const double& m, const double& Te)
// n: plasma density (m^-3)
// m: ion mass (kg)
// Te: electron temperature (K)
{
    return sqrt(m*pow(sig,2.)/(cts::kB*Te));
}

// lande g factor
double getLandeGFac(const double j, const double l, const double s)
// j: total angular momentum quantum number: J = L+S
// l: orbital angular momentum quantum number: L, note l>0
// s: spin angular momentum quantum number: S
{ return 1.+(j*(j+1.)+s*(s+1.)-l*(l+1.))/(2.*j*(j+1.)); }

// returns Zeeman shift in SI units
double getZeemanShift(const double B,const double gL, const double m)
// B: magnitude of magnetic field in Tesla
// gL: Lande g-factor (dimensionless)
// m: projection quantum number for total angular momentum J=L+S
{ return gL*cts::uB*B*m; }