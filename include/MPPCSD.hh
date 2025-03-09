#ifndef MPPCSD_HH
#define MPPCSD_HH

#include "G4VSensitiveDetector.hh"
#include "MPPCHit.hh"

#include "TSpline.h"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;

class MPPCSD : public G4VSensitiveDetector {
public:
  MPPCSD(const G4String& name);
  ~MPPCSD() override;

  void Initialize(G4HCofThisEvent* HCE) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  void EndOfEvent(G4HCofThisEvent* HCE) override;

private:
  G4THitsCollection<MPPCHit>* m_hits_collection;
  TSpline3* m_qe_spline;
  G4double m_range_min;
  G4double m_range_max;

  void InitializeQESpline();
};

#endif
