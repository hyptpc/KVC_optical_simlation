#include "EventAction.hh"
#include "AnaManager.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
}

EventAction::EventAction() {
}

EventAction::~EventAction() {
}

void EventAction::BeginOfEventAction(const G4Event* anEvent) {
  gAnaMan.BeginOfEventAction(anEvent);
}

void EventAction::EndOfEventAction(const G4Event* anEvent) {
  G4int eventID = anEvent->GetEventID();
  gAnaMan.EndOfEventAction(anEvent);

  // if (eventID % 100 == 0) {
    G4cout << "   Event number = " << eventID << G4endl;
  // }
}
