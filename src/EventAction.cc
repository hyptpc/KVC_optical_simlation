#include "EventAction.hh"
#include "AnaManager.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
}

EventAction::EventAction()
  : G4UserEventAction(),
    fNCherenkovGen(0),
    fNDeltaElectrons(0) {
}

EventAction::~EventAction() {
}

//_____________________________________________________________________________
void EventAction::BeginOfEventAction(const G4Event* anEvent) {
  gAnaMan.BeginOfEventAction(anEvent);
  fNCherenkovGen = 0; // Initialize Cherenkov photon count
  fNDeltaElectrons = 0; // Initialize delta electron count
}

//_____________________________________________________________________________
void EventAction::EndOfEventAction(const G4Event* anEvent) {
  G4int eventID = anEvent->GetEventID();

  gAnaMan.SetCherenkovGen(fNCherenkovGen);
  gAnaMan.SetNumDeltaElectrons(fNDeltaElectrons);
  gAnaMan.EndOfEventAction(anEvent); // Save event data to AnaManager


  if (eventID % 100 == 0) {
    G4cout << "   Event number = " << eventID << G4endl;
  }


}