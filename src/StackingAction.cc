// -*- C++ -*-

#include "StackingAction.hh"

#include <G4Electron.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>
#include <G4ios.hh>

#include "AnaManager.hh"
#include "EventAction.hh"
#include "KVC_TrackInfo.hh"

#define DEBUG 0

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();

  //___________________________________________________________________________
  EventAction*
  GetEventAction()
  {
    return static_cast<EventAction*>(
      G4EventManager::GetEventManager()->GetUserEventAction());
  }

  //___________________________________________________________________________
  G4bool
  IsInQuartz(const G4Track* track)
  {
    const G4VPhysicalVolume* volume = track->GetVolume();
    return (volume && volume->GetName() == "KvcPV");
  }
}

//_____________________________________________________________________________
StackingAction::StackingAction()
  : G4UserStackingAction(),
    m_n_scintillation_all(0),
    m_n_cerenkov_all(0),
    m_n_cerenkov_quartz(0)
{
}

//_____________________________________________________________________________
StackingAction::~StackingAction()
{
}

//_____________________________________________________________________________
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  // --- Optical photons ---
  if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() &&
      aTrack->GetParentID() > 0) { // secondary photon
    const auto creator = aTrack->GetCreatorProcess();
    if (!creator) return fUrgent;

    const G4String& process_name = creator->GetProcessName();
    if (process_name == "Scintillation") {
      ++m_n_scintillation_all;
    } else if (process_name == "Cerenkov") {
      ++m_n_cerenkov_all;

#if DEBUG
      const G4VPhysicalVolume* volume = aTrack->GetVolume();
      if (volume) {
        G4cout << "Cerenkov photon generated in volume: "
               << volume->GetName() << G4endl;
      } else {
        G4cout << "Cerenkov photon generated in an unknown volume" << G4endl;
      }
#endif

      const G4bool is_in_quartz = IsInQuartz(aTrack);
      const G4double energy = aTrack->GetKineticEnergy();
      if (is_in_quartz) {
        ++m_n_cerenkov_quartz;
        gAnaMan.AddGenWavelength((CLHEP::h_Planck * CLHEP::c_light / energy) / CLHEP::nm);
      }

      // Energy range of photons counted as generated Cherenkov photons
      constexpr G4double energy_min = 1.37 * eV;
      constexpr G4double energy_max = 3.87 * eV;

      if (is_in_quartz && energy >= energy_min && energy < energy_max) {
        auto event_action = GetEventAction();
        if (event_action) event_action->AddCherenkovGen();

        // Tag this track as "From Quartz"
        aTrack->SetUserInformation(new KVC_TrackInfo(true));
      }
    }
  }

  // --- Delta electrons (secondary electrons from ionization) ---
  // Count only new tracks: an electron suspended after emitting Cherenkov photons
  // (SetCerenkovTrackSecondariesFirst) comes back here every time it is resumed.
  if (aTrack->GetDefinition() == G4Electron::ElectronDefinition() &&
      aTrack->GetParentID() > 0 && aTrack->GetCurrentStepNumber() == 0) {
    const auto creator = aTrack->GetCreatorProcess();
    if (creator) {
      const G4String& process_name = creator->GetProcessName();
      if (process_name == "eIoni" || process_name == "ionIoni" || process_name == "hIoni") {
        if (IsInQuartz(aTrack)) {
          auto event_action = GetEventAction();
          if (event_action) event_action->AddDeltaElectron();
        }
      }
    }
  }

  return fUrgent;
}

//_____________________________________________________________________________
void
StackingAction::NewStage()
{
  gAnaMan.SetNumOfCerenkovAll(m_n_cerenkov_all);
  gAnaMan.SetNumOfCerenkovQuartz(m_n_cerenkov_quartz);
}

//_____________________________________________________________________________
void
StackingAction::PrepareNewEvent()
{
  m_n_scintillation_all = 0;
  m_n_cerenkov_all      = 0;
  m_n_cerenkov_quartz   = 0;
}
