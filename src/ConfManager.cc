// -*- C++ -*-

#include "ConfManager.hh"

#include <exception>
#include <fstream>
#include <sstream>
#include <string>

#include <G4ios.hh>

//_____________________________________________________________________________
ConfManager&
ConfManager::GetInstance()
{
  static ConfManager s_instance;
  return s_instance;
}

//_____________________________________________________________________________
ConfManager::ConfManager()
{
}

//_____________________________________________________________________________
void
ConfManager::Set(const G4String& key, const G4String& value)
{
  m_config_map[key] = value;
}

//_____________________________________________________________________________
G4bool
ConfManager::Check(const G4String& key) const
{
  return m_config_map.find(key) != m_config_map.end();
}

//_____________________________________________________________________________
G4String
ConfManager::Get(const G4String& key) const
{
  auto itr = m_config_map.find(key);
  if (itr != m_config_map.end()) {
    return itr->second;
  }
  G4cerr << "Warning: Config key '" << key << "' not found!" << G4endl;
  return "";
}

//_____________________________________________________________________________
G4double
ConfManager::GetDouble(const G4String& key) const
{
  auto itr = m_config_map.find(key);
  if (itr == m_config_map.end()) {
    G4cerr << "Warning: Config key '" << key << "' not found! Using 0.0" << G4endl;
    return 0.0;
  }
  try {
    return std::stod(itr->second);
  } catch (const std::exception& e) {
    G4cerr << "Error: Config key '" << key << "' has invalid double value '"
           << itr->second << "': " << e.what() << G4endl;
    throw;
  }
}

//_____________________________________________________________________________
G4int
ConfManager::GetInt(const G4String& key) const
{
  auto itr = m_config_map.find(key);
  if (itr == m_config_map.end()) {
    G4cerr << "Warning: Config key '" << key << "' not found! Using 0" << G4endl;
    return 0;
  }
  try {
    return std::stoi(itr->second);
  } catch (const std::exception& e) {
    G4cerr << "Error: Config key '" << key << "' has invalid int value '"
           << itr->second << "': " << e.what() << G4endl;
    throw;
  }
}

//_____________________________________________________________________________
G4String
ConfManager::GetPath(const G4String& key) const
{
  const G4String path = Get(key);
  if (path.empty() || path[0] == '/') return path;
  return m_conf_dir + path;
}

//_____________________________________________________________________________
void
ConfManager::LoadConfigFile(const G4String& filename)
{
  std::ifstream file(filename);
  if (!file) {
    G4cerr << "Error: Cannot open config file " << filename << G4endl;
    return;
  }

  const auto slash_pos = filename.find_last_of('/');
  m_conf_dir = (slash_pos == std::string::npos) ? "" : filename.substr(0, slash_pos + 1);

  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string key, value;
    if (iss >> key >> value) {
      m_config_map[key] = value;
    }
  }
}
