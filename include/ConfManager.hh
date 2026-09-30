// -*- C++ -*-

#ifndef CONF_MANAGER_HH
#define CONF_MANAGER_HH

#include <map>

#include <G4String.hh>
#include <G4Types.hh>

//_____________________________________________________________________________
class ConfManager
{
public:
  static ConfManager& GetInstance();

  G4String Get(const G4String& key) const;
  G4double GetDouble(const G4String& key) const;
  G4int    GetInt(const G4String& key) const;
  // Returns the value as a path; a relative path is resolved against the
  // directory of the loaded config file. Absolute paths are returned as is.
  G4String GetPath(const G4String& key) const;

  void   Set(const G4String& key, const G4String& value);
  void   LoadConfigFile(const G4String& filename);
  G4bool Check(const G4String& key) const;

private:
  ConfManager();

private:
  std::map<G4String, G4String> m_config_map;
  G4String m_conf_dir; // Directory of the loaded config file (with trailing '/')
};

#endif
