// -*- C++ -*-

#include "EventAction.hh"

#include <G4Event.hh>
#include <G4ios.hh>

#include "AnaManager.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
}

//_____________________________________________________________________________
EventAction::EventAction()
  : G4UserEventAction(),
    m_n_cherenkov_gen(0),
    m_n_delta_electrons(0)
{
}

//_____________________________________________________________________________
EventAction::~EventAction()
{
}

//_____________________________________________________________________________
void
EventAction::BeginOfEventAction(const G4Event* anEvent)
{
  gAnaMan.BeginOfEventAction(anEvent);
  m_n_cherenkov_gen = 0;
  m_n_delta_electrons = 0;
}

//_____________________________________________________________________________
void
EventAction::EndOfEventAction(const G4Event* anEvent)
{
  const G4int event_id = anEvent->GetEventID();

  gAnaMan.SetCherenkovGen(m_n_cherenkov_gen);
  gAnaMan.SetNumDeltaElectrons(m_n_delta_electrons);
  gAnaMan.EndOfEventAction(anEvent); // Save event data to AnaManager

  if (event_id % 100 == 0) {
    G4cout << "   Event number = " << event_id << G4endl;
  }
}
