// -*- C++ -*-

#include "SteppingAction.hh"

#include <G4OpticalPhoton.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>

#include "AnaManager.hh"
#include "KVC_TrackInfo.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
}

//_____________________________________________________________________________
SteppingAction::SteppingAction()
  : G4UserSteppingAction(),
    m_air_pv(nullptr),
    m_wrap_pv(nullptr)
{
}

//_____________________________________________________________________________
SteppingAction::~SteppingAction()
{
}

//_____________________________________________________________________________
void
SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track = aStep->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;

  // Cache physical volume pointers once (pointer comparison is much faster than string comparison)
  if (!m_air_pv) {
    auto pv_store = G4PhysicalVolumeStore::GetInstance();
    m_air_pv  = pv_store->GetVolume("KvcMotherPV", false);
    m_wrap_pv = pv_store->GetVolume("WrapPV", false);
  }

  // Photon detection is handled in MPPCSD. Here we only monitor lost photons.
  // A photon is "trapped/lost" if it is killed in the air or wrapper volumes.
  // Only photons born in the quartz (tagged with KVC_TrackInfo) are counted.
  if (track->GetTrackStatus() == fStopAndKill) {
    auto info = static_cast<KVC_TrackInfo*>(track->GetUserInformation());
    if (info && info->IsFromQuartz()) {
      const G4VPhysicalVolume* pre_pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
      if (pre_pv == m_air_pv || pre_pv == m_wrap_pv) {
        gAnaMan.IncrementTrappedAir();
      }
    }
  }
}
