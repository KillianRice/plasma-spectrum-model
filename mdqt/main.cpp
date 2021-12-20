// Cluster modules to load:
// module load iccifort/2020.1.217 impi/2019.7.217 Armadillo/9.880.1

// Compile and run on cluster:
// g++ -I source/header/ -std=c++2a -fopenmp main.cpp MDQT.cpp source/md/* source/qt/* source/wigner-symbols/* -o main -03 -lm -larmadillo -lstdc++fs

#include "PlasmaSettings.hpp"
#include <filesystem>

// forward declaration for MDQT program
void MDQT(PlasmaSettings pms,int array,std::filesystem::path savePath);

int main(int argc, char *argv[]) 
{
    // process command line arguments
    std::filesystem::path save_path {getCommandLineArg(argc,argv, "-p","--path")}; // relative top-level save path for data
    std::filesystem::create_directories(save_path); 

    std::filesystem::path settings_path {getCommandLineArg(argc,argv, "-s","--settings")}; // relative path to settings file

    int task_array;
    std::string task_array_str = getCommandLineArg(argc, argv, "-a", "--array"); // read in task array
    if(task_array_str.empty()) task_array = 0; // default array index is zero when array not specified
    else task_array = std::stoi(task_array_str); // use input task array if provided

    // load .settings file
    PlasmaSettings pms(settings_path,task_array);

    // run MDQT program
    std::cout << "Task array range: [0 " << pms.mat().size()-1 << "]" << std::endl;
    int set = floor(task_array/pms.runs());
    save_path /= "set" + num2str(set);
    std::filesystem::create_directories(save_path);
    save_path /= "run" + num2str(pms.getvar("runs"));
    std::filesystem::create_directories(save_path);
    pms.write_array_params(save_path);
    MDQT(pms,task_array,save_path);

    return 0;
}