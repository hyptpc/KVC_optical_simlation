#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class EventAction : public G4UserEventAction {

    //追加
private: 
    G4int fNCerGen;  // チェレンコフ光生成数
    //追加終

public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    //追加
    void AddCerGen() { fNCerGen++; } //addCerGenは、fNCerGen++という操作
    G4int GetCerGen() const { return fNCerGen; }
    //追加終
    
};

#endif
