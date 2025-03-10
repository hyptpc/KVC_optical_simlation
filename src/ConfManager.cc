#include "ConfManager.hh"
#include <fstream>
#include <sstream>
#include <iostream>

ConfManager& ConfManager::Instance() {
    static ConfManager instance;
    return instance;
}

ConfManager::ConfManager() {}

void ConfManager::Set(const std::string& key, const std::string& value) {
    config_map[key] = value;
}

std::string ConfManager::Get(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return it->second;
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return "";
}

double ConfManager::GetDouble(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return std::stod(it->second);
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return 0.0;
}

int ConfManager::GetInt(const std::string& key) const {
    auto it = config_map.find(key);
    if (it != config_map.end()) {
        return std::stoi(it->second);
    }
    std::cerr << "Warning: Config key '" << key << "' not found!" << std::endl;
    return 0;
}

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
