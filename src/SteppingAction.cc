// -*- C++ -*-

#include "SteppingAction.hh"

#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
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
    m_kvc_pv(nullptr),
    m_mppc_lv(nullptr)
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

  // Cache volume pointers once (pointer comparison is much faster than string comparison)
  if (!m_kvc_pv) {
    m_kvc_pv  = G4PhysicalVolumeStore::GetInstance()->GetVolume("KvcPV", false);
    m_mppc_lv = G4LogicalVolumeStore::GetInstance()->GetVolume("MppcLV", false);
  }

  // Photon detection is handled in MPPCSD. Here we only monitor lost photons.
  // A photon born in the quartz (tagged with KVC_TrackInfo) is "lost" if it is
  // killed outside the quartz (absorbed in the air, wrapper, blacksheet, ... or
  // leaving the world), excluding photons reaching the MPPCs.
  if (track->GetTrackStatus() == fStopAndKill) {
    auto info = static_cast<KVC_TrackInfo*>(track->GetUserInformation());
    if (info && info->IsFromQuartz()) {
      const G4VPhysicalVolume* pre_pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
      const G4bool is_in_quartz = (pre_pv == m_kvc_pv);
      const G4bool is_in_mppc   = (pre_pv && pre_pv->GetLogicalVolume() == m_mppc_lv);
      if (!is_in_quartz && !is_in_mppc) {
        gAnaMan.IncrementTrappedAir();
      }
    }
  }
}
