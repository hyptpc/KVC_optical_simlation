#ifndef CONFMANAGER_HH
#define CONFMANAGER_HH

#include <string>
#include <unordered_map>

class ConfManager {
public:
    static ConfManager& Instance();
    
    std::string Get(const std::string& key) const;
    double GetDouble(const std::string& key) const;
    int GetInt(const std::string& key) const;

    void Set(const std::string& key, const std::string& value);
    void LoadConfigFile(const std::string& filename);

private:
    ConfManager();
    std::unordered_map<std::string, std::string> config_map;
};

#endif // CONFMANAGER_HH
