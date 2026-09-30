// -*- C++ -*-

#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

//_____________________________________________________________________________
ActionInitialization::ActionInitialization()
  : G4VUserActionInitialization()
{
}

//_____________________________________________________________________________
ActionInitialization::~ActionInitialization()
{
}

//_____________________________________________________________________________
void
ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

//_____________________________________________________________________________
void
ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction);
  SetUserAction(new EventAction);
  SetUserAction(new SteppingAction);
  SetUserAction(new StackingAction);
}
