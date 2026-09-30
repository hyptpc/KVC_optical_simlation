#ifndef CONFMANAGER_HH
#define CONFMANAGER_HH

#include <string>
#include <unordered_map>

class ConfManager {
public:
    static ConfManager& GetInstance();
    
    std::string Get(const std::string& key) const;
    double GetDouble(const std::string& key) const;
    int GetInt(const std::string& key) const;
    // Returns the value as a path; a relative path is resolved against the
    // directory of the loaded config file. Absolute paths are returned as is.
    std::string GetPath(const std::string& key) const;

    void Set(const std::string& key, const std::string& value);
    void LoadConfigFile(const std::string& filename);
    bool Check(const std::string& key) const;

private:
    ConfManager();
    std::unordered_map<std::string, std::string> config_map;
    std::string conf_dir; // Directory of the loaded config file (with trailing '/')
};

#endif // CONFMANAGER_HH
