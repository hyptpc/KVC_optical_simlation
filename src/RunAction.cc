#include "RunAction.hh"
#include "AnaManager.hh"

#include <fstream>

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4StateManager.hh>
#include <G4Timer.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "AnaManager.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
G4Timer timer;
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
  G4Random::setTheSeed(std::time(nullptr));
  timer.Start();
}

//_____________________________________________________________________________
void
RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer.Stop();
  gAnaMan.EndOfRunAction(aRun);
  G4cout << "   Process end  = " << timer.GetClockTime()
	 << "   Event number = " << aRun->GetNumberOfEvent() << G4endl
	 << "   Elapsed time = " << timer << G4endl << G4endl;
}
