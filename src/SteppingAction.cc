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

SteppingAction::SteppingAction() {
}

SteppingAction::~SteppingAction()
{}

// void SteppingAction::UserSteppingAction(const G4Step* step)
// {
// }


void SteppingAction::UserSteppingAction(const G4Step* step)
{
}
