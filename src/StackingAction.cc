#include "StackingAction.hh"

#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ClassificationOfNewTrack.hh"

#include "AnaManager.hh"

#include "G4EventManager.hh"
#include "EventAction.hh"

#include "G4SystemOfUnits.hh"      
#include "G4PhysicalConstants.hh"  

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
}

#define DEBUG 0

//_____________________________________________________________________________
StackingAction::StackingAction()
  : G4UserStackingAction(),
    fScintillationAll(0), fCerenkovAll(0), fCerenkovQuartz(0)
{}

//_____________________________________________________________________________
StackingAction::~StackingAction()
{}

//_____________________________________________________________________________
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{    

  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) // Focus on optical photons
  { 
    if(aTrack->GetParentID() > 0){ // particle is secondary
      const auto* creator = aTrack->GetCreatorProcess();
      if (!creator) return fUrgent;     

      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation") // Scintillation photon
        ++fScintillationAll;
      else if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov") { // Cerenkov photon
	      ++fCerenkovAll;

	  const G4VPhysicalVolume* volume = aTrack->GetVolume(); // Current volume of the photon
    const bool in_quartz = (volume && volume->GetName() == "KvcPV");
    if (in_quartz) ++fCerenkovQuartz;

    constexpr G4double Emin = 1.37 * eV;
    constexpr G4double Emax = 3.87 * eV;
    const G4double E = aTrack->GetKineticEnergy();

    if(in_quartz && E >= Emin && E < Emax ){
    auto eventAction = static_cast<EventAction*>(
    G4EventManager::GetEventManager()->GetUserEventAction());
    if (eventAction) eventAction->AddCherenkovGen(); // Increment Cherenkov count
    }
    
  }

#if DEBUG
	if (volume) {
	  G4cout << "Cerenkov photon generated in volume: " 
		 << volume->GetName() << G4endl;
	} else {
	  G4cout << "Cerenkov photon generated in an unknown volume" << G4endl;
	}
#endif
        
      }
    }
  
	
  return fUrgent;
}


//_____________________________________________________________________________
void StackingAction::NewStage()
{
  // G4cout << "Number of Scintillation photons produced in this event : "
  // 	 << fScintillationAll << G4endl;
  // G4cout << "Number of Cerenkov photons produced in this event : "
  // 	 << fCerenkovAll << G4endl;
  // G4cout << "Number of Cerenkov photons produced in Quartz : "
  // 	 << fCerenkovQuartz << G4endl;
  gAnaMan.SetNumOfCerenkovAll(fCerenkovAll);
  gAnaMan.SetNumOfCerenkovQuartz(fCerenkovQuartz);
}

//_____________________________________________________________________________
void StackingAction::PrepareNewEvent()
{
  fScintillationAll = 0;
  fCerenkovAll      = 0;
  fCerenkovQuartz   = 0;
}
