#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class EventAction : public G4UserEventAction {

    //追加
private: 
    G4int fNCherenkovGen;  // チェレンコフ光生成数
    //追加終

public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    //追加
    void AddCherenkovGen() { fNCherenkovGen++; } //addCherenkovGenは、fNCherenkovGen++という操作
    G4int GetCherenkovGen() const { return fNCherenkovGen; }
    //追加終
    
};

#endif
