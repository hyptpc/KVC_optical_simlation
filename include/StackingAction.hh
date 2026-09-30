// -*- C++ -*-

#ifndef STACKING_ACTION_HH
#define STACKING_ACTION_HH

#include <G4Types.hh>
#include <G4UserStackingAction.hh>

class G4Track;

//_____________________________________________________________________________
class StackingAction : public G4UserStackingAction
{
public:
  StackingAction();
  ~StackingAction() override;

public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) override;
  void NewStage() override;
  void PrepareNewEvent() override;

private:
  G4int m_n_scintillation_all; // Scintillation photons in this event (not stored yet)
  G4int m_n_cerenkov_all;      // Cherenkov photons in this event (all volumes)
  G4int m_n_cerenkov_quartz;   // Cherenkov photons generated in the quartz
};

#endif
