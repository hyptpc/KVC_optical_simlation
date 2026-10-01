// -*- C++ -*-

#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>

class G4LogicalVolume;
class G4Step;
class G4VPhysicalVolume;

//_____________________________________________________________________________
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  ~SteppingAction() override;

public:
  void UserSteppingAction(const G4Step* aStep) override;

private:
  G4VPhysicalVolume* m_kvc_pv;  // Quartz radiator
  G4LogicalVolume*   m_mppc_lv; // MPPCs
};

#endif
