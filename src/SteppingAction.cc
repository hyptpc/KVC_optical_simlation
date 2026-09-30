#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "AnaManager.hh"
#include "KVC_TrackInfo.hh"

SteppingAction::SteppingAction()
  : fAirVol(nullptr), fWrapVol(nullptr)
{
}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  if(track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;

  // Cache physical volume pointers once (Pointer comparison is MUCH faster than string comparison)
  if(!fAirVol) {
      auto pvStore = G4PhysicalVolumeStore::GetInstance();
      fAirVol  = pvStore->GetVolume("KvcMotherPV", false);
      fWrapVol = pvStore->GetVolume("WrapPV", false);
  }

  // Photon detection is handled in MPPCSD. Here we only monitor lost photons.

  // --- Monitoring Logics ---
  // A photon is "trapped/lost" if it is KILLED in the Air or Wrap volumes.
  // We ONLY count photons that were born in Quartz (checked via TrackInfo).
  if (track->GetTrackStatus() == fStopAndKill) {
      auto info = static_cast<KVC_TrackInfo*>(track->GetUserInformation());
      if (info && info->IsFromQuartz()) {
          G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
          if(preVol == fAirVol || preVol == fWrapVol) {
              AnaManager::GetInstance().IncrementTrappedAir();
          }
      }
  }
}
