// -*- C++ -*-

#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include <G4UserRunAction.hh>

class G4Run;

//_____________________________________________________________________________
class RunAction : public G4UserRunAction
{
public:
  RunAction();
  ~RunAction() override;

public:
  void BeginOfRunAction(const G4Run* aRun) override;
  void EndOfRunAction(const G4Run* aRun) override;
};

#endif
