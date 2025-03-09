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

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
}

StackingAction::StackingAction()
  : G4UserStackingAction(),
    fScintillationAll(0), fCerenkovAll(0), fCerenkovQuartz(0)
{}

StackingAction::~StackingAction()
{}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{    

  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  {  // particle is optical photon
    if(aTrack->GetParentID() > 0)
    {  // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        ++fScintillationAll;
      else if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov") {
	++fCerenkovAll;

	const G4VPhysicalVolume* volume = aTrack->GetVolume();
	if (volume->GetName() == "KvcPV") ++fCerenkovQuartz;
	
	// if (volume) {
	//   G4cout << "Cerenkov photon generated in volume: " 
	// 	 << volume->GetName() << G4endl;
	// } else {
	//   G4cout << "Cerenkov photon generated in an unknown volume" << G4endl;
	// }
        
      }
    }
  }
	
  return fUrgent;
}


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

void StackingAction::PrepareNewEvent()
{
  fScintillationAll = 0;
  fCerenkovAll      = 0;
  fCerenkovQuartz   = 0;
}
