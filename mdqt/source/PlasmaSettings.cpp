#include "PlasmaSettings.hpp"

// class constructor
PlasmaSettings::PlasmaSettings(std::filesystem::path settings_path, int array):
    m_array{array}
{
    load_settings(settings_path);
    check_units();
    process_runs();

    if (array < 0 || array > m_mat.size()) std::cerr << "User-specified <array> value must correspond to row index of m_unique_val.";
    m_val = m_mat[m_array];
}

// load .settings file that defines plasma characteristics for MHD simulation
void PlasmaSettings::load_settings(const std::filesystem::path& settings_path)
{
    // mat_in[i] contains vector with [name unit values]
    std::vector<std::vector<std::string>> mat_in = readCSV(settings_path);
    std::vector<std::vector<std::string>> vars(mat_in.size());
    for (int i = 0; i < mat_in.size(); i++){
        for (int j = 2; j < mat_in[i].size(); j++){
            vars[i].push_back(mat_in[i][j]);
        }
    }

    // get unique combinations of vars -> store in <pms>
    m_names.resize(mat_in.size());
    m_units.resize(mat_in.size());
    for (int i = 0; i < mat_in.size(); i++) m_names[i] = mat_in[i][0];
    for (int i = 0; i < mat_in.size(); i++) m_units[i] = mat_in[i][1];
    m_mat = unique_comb(vars[0],vars[1]);
    for (int i = 2; i < vars.size(); i++) m_mat = unique_comb(m_mat,vars[i]);
}

// verify that all variable units are specified correctly: either <SI>, <opt>, or another variable name
void PlasmaSettings::check_units(void) const{
    // check that units are same size as varnames
    bool size_check = m_names.size() == m_units.size();
    if (!size_check) std::cerr << "<varnames> and <varunits> must be initialized with the same size";

    // possible_units are "SI" or a variable name
    std::vector<std::string> possible_units;
    possible_units.reserve(m_names.size()+2);
    possible_units.push_back("SI");
    possible_units.push_back("opt");
    for (auto& i : m_names) possible_units.push_back(i);

    // check that specified units are possible
    for (auto& i : m_units){
        bool unit_check{false};
        for (auto& j: possible_units){
            if (j.find(i) == 0) unit_check = true;
        }
        if (!unit_check) std::cerr << "Specified unit <" <<  i << "> must be <SI> or another variable name." << std::endl;
    } 

    // for varialbes with non-SI units, check that their dependencies are in SI
    for (int i = 0; i < m_units.size(); i++){ // for each unit
        if (!strcmp(m_units[i],"SI") && !strcmp(m_units[i],"opt")){ // if unit is not "SI"
            // find which variable m_units[i] is in terms of 
            int dependent;
            for (int j = 0; j < m_names.size(); j++){
                if (strcmp(m_names[j],m_units[i])){
                    dependent = j;
                    break;
                }
            }
            // if the dependent is not in "SI", throw an error
            if (!strcmp(m_units[dependent],"SI")) std::cerr << "Variable <" << m_names[i] << "> is specified in terms of <" << m_names[dependent] << ">, which is not specified in <SI> units." << std::endl;
        }
    }
    
}

// return numeric variable with SI units
double PlasmaSettings::getvar(const std::string& name) const
{
    // ensure given name corresponds to variable name and that variable is not type "opt"
    bool found{false};
    int iter;
    for (int i{0}; i < m_names.size(); i++){ 
        if (strcmp(m_names[i],name)){ 
            found = true;
            iter = i;
            break;
        }
    }
    if (!found) std::cerr << "User input <" << name << "> does not match a variable name." << std::endl;
    if (strcmp(m_units[iter],"opt")) std::cerr << "Requested variable <" << name << "> is of type <opt>." << std::endl;

    // get value of variable corresponding to <name>
    double val = stod(m_val[iter]);

    // if not in SI units, convert <val> to SI units
    if (!strcmp(m_units[iter],"SI")){ // if variable units are not "SI" (i.e., expressed in terms of another variable)
        for (int i = 0; i < m_names.size(); i ++){ // for each variable
            if (strcmp(m_names[i],m_units[iter])){ // if name matches m_units[iter]
                val *= stod(m_val[i]); // record value as conversion factor
            }
        } 
    }

    return val;
}

// return option variable as a string
std::string PlasmaSettings::getopt(const std::string& name) const
{
    // ensure given name has units "opt"
    bool found{false};
    int iter;
    for (int i{0}; i < m_names.size(); i++){
        if (strcmp(m_names[i],name)){
            found = true;
            iter = i;
            break;
        }
    }
    if (!found) std::cerr << "User input <" << name << "> does not match a variable name." << std::endl;
    if (!strcmp(m_units[iter],"opt")) std::cerr << "Requested variable <" << name << "> is not of type <opt>." << std::endl;
    return m_val[iter];



}

// write plasma settings for specified m_array value (i.e., the corresponding row of m_mat)
void PlasmaSettings::write_array_params(const std::filesystem::path& path) const
{
    std::filesystem::create_directory(path);
    std::filesystem::path file_path{path/"plasma.settings"};
    if (std::filesystem::exists(file_path)) std::ofstream out_file(file_path);
    append_row_to_csv(file_path,m_names);
    append_row_to_csv(file_path,m_units);
    append_row_to_csv(file_path,m_val);
}

std::vector<std::vector<std::string>> PlasmaSettings::mat() const{ return m_mat; }
int PlasmaSettings::runs() const {return m_runs;}

// create unique run numbers
void PlasmaSettings::process_runs()
{
    bool found = false;
    size_t loc;
    for (int i = 0; i < m_names.size(); i++){
        if (strcmp(m_names[i],"runs")){
            found = true;
            loc = i;
            break;
        } 
    }
    if (!found) std::cerr << ".settings file must have <runs> specified." << std::endl;

    std::vector<std::vector<std::string>> new_mat;
    m_runs = stod(m_mat[0][loc]);
    new_mat.reserve(m_mat.size()*m_runs);
    for (int i = 0; i < m_mat.size(); i++){
        for (int j = 0; j < m_runs; j++){
            new_mat.push_back(m_mat[i]);
            new_mat.back()[loc] = num2str(j);
        }
    }
    m_mat = new_mat;
}