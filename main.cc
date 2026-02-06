#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "AnaManager.hh"
#include "RunAction.hh"
#include "ConfManager.hh"
    
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManager.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4Cerenkov.hh"
#include "G4DecayPhysics.hh"

#include <random>

namespace
{
  auto& gAnaMan  = AnaManager::GetInstance();
  auto& gConfMan = ConfManager::GetInstance();
  void PrintUsage()
  {
    G4cerr << " Usage: " << G4endl
	   << " KVCOpticalSim <conf file> <output rootfile name> [macro]"
           << G4endl;
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc < 3 || argc > 4) {
    PrintUsage();
    return 1;
  }
  gConfMan.LoadConfigFile(argv[1]); 
  gAnaMan.SetOutputRootfilePath(argv[2]);
  
  G4String macro;
  if (argc == 4) macro = argv[3];

  G4UIExecutive* ui = nullptr;
  if (macro.empty())
  {
    ui = new G4UIExecutive(argc, argv);
  }

  auto runManager = new G4RunManager();

  std::random_device rd;
  long seed;
  if (gConfMan.Check("seed")) {
      seed = gConfMan.GetInt("seed");
      G4cout << "Random seed: " << seed << " (Fixed from config)" << G4endl;
  } else {
      seed = rd();
      G4cout << "Random seed: " << seed << " (Randomized)" << G4endl;
  }
  G4Random::setTheSeed(seed);

  runManager->SetUserInitialization(new DetectorConstruction());

  // Physics List setting
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  auto opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  if (gConfMan.GetInt("decay") == 1) physicsList->RegisterPhysics(new G4DecayPhysics());
  runManager->SetUserInitialization(physicsList);

  // G4Cerenkov setting
  auto optical_params = G4OpticalParameters::Instance();
  optical_params->SetCerenkovMaxPhotonsPerStep(100);
  optical_params->SetCerenkovStackPhotons(true);
  optical_params->SetCerenkovTrackSecondariesFirst(true);
  optical_params->SetCerenkovVerboseLevel(1);
  optical_params->SetBoundaryVerboseLevel(1);
  optical_params->SetAbsorptionVerboseLevel(1);
    
  runManager->SetUserInitialization(new ActionInitialization());
  runManager->Initialize();

  
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();


  if (!macro.empty())
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}





