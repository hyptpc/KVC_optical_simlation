#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "AnaManager.hh"
#include "RunAction.hh"

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

namespace
{
  void PrintUsage()
  {
    G4cerr << " Usage: " << G4endl;
#ifdef GEANT4_USE_GDML
    G4cerr << " MPPCSim [-g gdmlfile] [-m macro ] [-u UIsession] [-r seed] "
           << G4endl;
#else
    G4cerr << " MPPCSim  [-m macro ] [-u UIsession] [-r seed] "
           << G4endl;
#endif
  }
}  // namespace

int main(int argc, char** argv)
{
  if (argc > 7)
  {
    PrintUsage();
    return 1;
  }

  G4String gdmlfile;
  G4String macro;
  G4String session;
  G4long myseed = 345354;

  for (G4int i = 1; i < argc; i = i + 2)
  {
    if (G4String(argv[i]) == "-g")
      gdmlfile = argv[i + 1];
    else if (G4String(argv[i]) == "-m")
      macro = argv[i + 1];
    else if (G4String(argv[i]) == "-u")
      session = argv[i + 1];
    else if (G4String(argv[i]) == "-r")
      myseed = atoi(argv[i + 1]);
    else
    {
      PrintUsage();
      return 1;
    }
  }

  G4UIExecutive* ui = nullptr;
  if (macro.empty())
  {
    ui = new G4UIExecutive(argc, argv);
  }

  auto runManager = new G4RunManager();

  // AnaManagerの初期化
  AnaManager::GetInstance();

  // 乱数シードの設定
  G4Random::setTheSeed(myseed);

  // Detector 初期化
  runManager->SetUserInitialization(new DetectorConstruction());

  // 物理リスト設定
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  auto opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  // physicsList->RegisterPhysics(new G4DecayPhysics());  
  runManager->SetUserInitialization(physicsList);

    
  auto optical_params = G4OpticalParameters::Instance();

  // G4Cerenkovの設定
  optical_params->SetCerenkovMaxPhotonsPerStep(100);
  optical_params->SetCerenkovStackPhotons(true);
  optical_params->SetCerenkovTrackSecondariesFirst(true);
  optical_params->SetCerenkovVerboseLevel(1);
  optical_params->SetBoundaryVerboseLevel(1);
  
  // // // G4Scintillationの設定
  // // optical_params->SetScintByParticleType(false);
  // // optical_params->SetScintTrackInfo(false);
  // // optical_params->SetScintTrackSecondariesFirst(true);
  // // optical_params->SetScintFiniteRiseTime(false);
  // // optical_params->SetScintStackPhotons(true);
  // // optical_params->SetScintVerboseLevel(1);

  // G4OpAbsorptionの設定
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
