// -*- C++ -*-

#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include <G4Types.hh>
#include <G4UserEventAction.hh>

class G4Event;

//_____________________________________________________________________________
class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

private:
  G4int m_n_cherenkov_gen;   // Number of generated Cherenkov photons
  G4int m_n_delta_electrons; // Number of generated delta electrons

public:
  void BeginOfEventAction(const G4Event* anEvent) override;
  void EndOfEventAction(const G4Event* anEvent) override;

  void  AddCherenkovGen() { ++m_n_cherenkov_gen; }
  G4int GetCherenkovGen() const { return m_n_cherenkov_gen; }
  void  AddDeltaElectron() { ++m_n_delta_electrons; }
  G4int GetDeltaElectrons() const { return m_n_delta_electrons; }
};

#endif
