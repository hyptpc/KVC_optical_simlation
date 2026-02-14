#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4TouchableHandle.hh"
#include "MPPCSD.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "KVC_OpticalProperties.hh"
#include "AnaManager.hh"
#include "KVC_TrackInfo.hh"
#include "ConfManager.hh"
#include "G4EventManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"

#define KVC_DEBUG_STEPPING 0

SteppingAction::SteppingAction() 
  : fOpProcess(nullptr), fAirVol(nullptr), fWrapVol(nullptr) 
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

  // Retrieve OpBoundaryProcess if not cached
  if(!fOpProcess){
    G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i=0; i<nprocesses; ++i){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        fOpProcess = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  if(!fOpProcess) return;

  G4OpBoundaryProcessStatus status = fOpProcess->GetStatus();
  G4bool detected = false;

#ifdef USE_SURFACE_PDE
  // --- Detection Logic (PDE check) ---
  if (status == Detection) {
      detected = true;
  } else if (status == FresnelRefraction) {
      // Check if it entered MPPC
      G4TouchableHandle touchable = step->GetPostStepPoint()->GetTouchableHandle();
      auto vol = touchable->GetVolume();
      if(vol && G4StrUtil::contains(vol->GetName(), "Mppc")){
          const G4LogicalSkinSurface* surf = G4LogicalSkinSurface::GetSurface(vol->GetLogicalVolume());
          if(surf) {
            G4SurfaceProperty* prop = const_cast<G4SurfaceProperty*>(surf->GetSurfaceProperty());
            G4OpticalSurface* optSurf = dynamic_cast<G4OpticalSurface*>(prop);
            if(optSurf){
                G4MaterialPropertiesTable* MPT = optSurf->GetMaterialPropertiesTable();
                if(MPT && MPT->GetProperty("EFFICIENCY")){
                    G4double energy = track->GetTotalEnergy();
                    G4double eff = MPT->GetProperty("EFFICIENCY")->Value(energy); 
                    G4double qe_scale = ConfManager::GetInstance().GetDouble("qe_scale");
                    G4double prob = eff * qe_scale;
                    if(prob > 1.0) prob = 1.0; 
                    if(G4UniformRand() < prob) {
                        detected = true;
                    }
                }
            }
          }
      }
  }
#endif

  if(detected){
    // Determine which volume detected it.
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    G4TouchableHandle touchable = postStepPoint->GetTouchableHandle();
    
    G4ThreeVector worldPos = postStepPoint->GetPosition();
    G4ThreeVector pos      = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    G4double hitTime = postStepPoint->GetGlobalTime();
    G4double energy = track->GetTotalEnergy();
    G4double waveLength = (CLHEP::h_Planck * CLHEP::c_light / energy) / CLHEP::nm; 
    G4int copyNumber = touchable->GetCopyNumber();
    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    MPPCHit* aHit = new MPPCHit();
    aHit->SetPosition(pos);
    aHit->SetWorldPosition(worldPos);
    aHit->SetEnergy(energy);
    aHit->SetWaveLength(waveLength);
    aHit->SetTime(hitTime);
    aHit->SetParticleID(track->GetDefinition()->GetPDGEncoding());
    aHit->SetCopyNumber(copyNumber);
    aHit->SetEventID(eventID);
    aHit->SetDetectFlag(1); // Detected!

    // Add to Collection
    auto HCTE = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetHCofThisEvent();
    if(HCTE){
       static G4int collID = -1;
       if(collID < 0) {
         collID = G4SDManager::GetSDMpointer()->GetCollectionID("MppcCollection");
       }
       if(collID >= 0){
          auto hitsCollection = (G4THitsCollection<MPPCHit>*)(HCTE->GetHC(collID));
          if(hitsCollection) hitsCollection->insert(aHit);
          else delete aHit;
       } else delete aHit;
    } else delete aHit;

    // IMPORTANT: Kill the track after detection!
    track->SetTrackStatus(fStopAndKill);
  }

  // --- Monitoring Logics ---
  // A photon is "trapped/lost" if it is KILLED in the Air or Wrap volumes, but NOT detected.
  // We ONLY count photons that were born in Quartz (checked via TrackInfo).
  if (!detected && track->GetTrackStatus() == fStopAndKill) {
      auto info = static_cast<KVC_TrackInfo*>(track->GetUserInformation());
      if (info && info->IsFromQuartz()) {
          G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
          if(preVol == fAirVol || preVol == fWrapVol) {
              AnaManager::GetInstance().IncrementTrappedAir();
          }
      }
  }
}
