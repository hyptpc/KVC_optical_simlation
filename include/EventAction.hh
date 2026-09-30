#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class EventAction : public G4UserEventAction {

private: 
    G4int fNCherenkovGen;    // Number of generated Cherenkov photons
    G4int fNDeltaElectrons;  // Number of generated delta electrons

public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddCherenkovGen() { fNCherenkovGen++; }
    G4int GetCherenkovGen() const { return fNCherenkovGen; }
    void AddDeltaElectron() { fNDeltaElectrons++; }
    G4int GetDeltaElectrons() const { return fNDeltaElectrons; }
    
};

#endif
