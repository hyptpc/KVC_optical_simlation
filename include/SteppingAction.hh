#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class G4Track;
class G4OpBoundaryProcess;
class G4VPhysicalVolume;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
  
private:
  G4OpBoundaryProcess* fOpProcess;
  G4VPhysicalVolume* fAirVol;
  G4VPhysicalVolume* fWrapVol;
};

#endif
