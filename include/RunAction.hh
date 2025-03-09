#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"

class RunAction : public G4UserRunAction {
public:
  RunAction();
  virtual ~RunAction();
  virtual void BeginOfRunAction(const G4Run* aRun);
  virtual void EndOfRunAction(const G4Run* aRun);
};

#endif
