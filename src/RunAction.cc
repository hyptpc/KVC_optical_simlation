// -*- C++ -*-

#include "RunAction.hh"

#include <G4Run.hh>
#include <G4Timer.hh>
#include <G4ios.hh>

#include "AnaManager.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
  G4Timer g_timer;
}

//_____________________________________________________________________________
RunAction::RunAction()
  : G4UserRunAction()
{
}

//_____________________________________________________________________________
RunAction::~RunAction()
{
}

//_____________________________________________________________________________
void
RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "   Run# = " << aRun->GetRunID() << G4endl;
  gAnaMan.BeginOfRunAction(aRun);
  g_timer.Start();
}

//_____________________________________________________________________________
void
RunAction::EndOfRunAction(const G4Run* aRun)
{
  g_timer.Stop();
  gAnaMan.EndOfRunAction(aRun);
  G4cout << "   Process end  = " << g_timer.GetClockTime()
         << "   Event number = " << aRun->GetNumberOfEvent() << G4endl
         << "   Elapsed time = " << g_timer << G4endl << G4endl;
}
