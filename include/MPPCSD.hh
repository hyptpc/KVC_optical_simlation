// -*- C++ -*-

#ifndef MPPC_SD_HH
#define MPPC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "MPPCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class TSpline3;

//_____________________________________________________________________________
class MPPCSD : public G4VSensitiveDetector
{
public:
  MPPCSD(const G4String& name);
  ~MPPCSD() override;

public:
  void   Initialize(G4HCofThisEvent* HCTE) override;
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) override;
  void   EndOfEvent(G4HCofThisEvent* HCTE) override;

private:
  void InitializeQESpline();

private:
  G4THitsCollection<MPPCHit>* m_hits_collection;
  TSpline3*                   m_qe_spline; // PDE as a function of photon energy
  G4double                    m_range_min; // Energy range of the PDE table
  G4double                    m_range_max;
  G4double                    m_qe_scale;  // Scale factor for the PDE (conf: qe_scale)
};

#endif
