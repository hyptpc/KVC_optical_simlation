#include "ConfManager.hh"
#include <fstream>
#include <sstream>
#include <iostream>

//_____________________________________________________________________________
ConfManager& ConfManager::GetInstance() {
    static ConfManager instance;
    return instance;
}

//_____________________________________________________________________________
ConfManager::ConfManager() {}

//_____________________________________________________________________________
void ConfManager::Set(const std::string& key, const std::string& value) {
    config_map[key] = value;
}

//_____________________________________________________________________________
bool ConfManager::Check(const std::string& key) const {
    return config_map.find(key) != config_map.end();
}

//_____________________________________________________________________________
std::string ConfManager::Get(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return it->second;
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return "";
}

//_____________________________________________________________________________
double ConfManager::GetDouble(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return std::stod(it->second);
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return 0.0;
}

//_____________________________________________________________________________
int ConfManager::GetInt(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return std::stoi(it->second);
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return 0;
}

//_____________________________________________________________________________
void ConfManager::LoadConfigFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open config file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (iss >> key >> value) {
            config_map[key] = value;
        }
    }
}
