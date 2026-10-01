// -*- C++ -*-

#include <random>

#include <G4BuilderType.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalParameters.hh>
#include <G4OpticalPhysics.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <QGSP_BERT.hh>
#include <Randomize.hh>

#include "ActionInitialization.hh"
#include "AnaManager.hh"
#include "ConfManager.hh"
#include "DetectorConstruction.hh"

namespace
{
  auto& gAnaMan  = AnaManager::GetInstance();
  auto& gConfMan = ConfManager::GetInstance();

  //___________________________________________________________________________
  void
  PrintUsage()
  {
    G4cerr << " Usage: " << G4endl
           << " KVCOpticalSim <conf file> <output rootfile name> [macro]"
           << G4endl;
  }
}

//_____________________________________________________________________________
int
main(int argc, char** argv)
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
  if (macro.empty()) {
    ui = new G4UIExecutive(argc, argv);
  }

  auto run_manager = new G4RunManager();

  // Random seed: fixed from the conf file if given, otherwise randomized
  std::random_device random_device;
  long seed;
  if (gConfMan.Check("seed")) {
    seed = gConfMan.GetInt("seed");
    G4cout << "Random seed: " << seed << " (Fixed from config)" << G4endl;
  } else {
    seed = random_device();
    G4cout << "Random seed: " << seed << " (Randomized)" << G4endl;
  }
  G4Random::setTheSeed(seed);

  run_manager->SetUserInitialization(new DetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physics_list = new QGSP_BERT;
  physics_list->ReplacePhysics(new G4EmStandardPhysics_option4());
  physics_list->RegisterPhysics(new G4OpticalPhysics());
  // Decay physics (G4DecayPhysics) is already included in QGSP_BERT.
  // Remove it when decay is disabled in the conf file (decay 0).
  if (gConfMan.GetInt("decay") == 0) physics_list->RemovePhysics(bDecay);
  // Production cut (range) for all particles; 0.1 mm if not given in the conf file.
  // 0.1 mm corresponds to an e- threshold of ~135 keV in quartz, below the
  // Cherenkov threshold of electrons (~0.19 MeV), so that delta rays emitting
  // Cherenkov light are produced.
  const G4double production_cut =
    (gConfMan.Check("production_cut") ? gConfMan.GetDouble("production_cut") : 0.1) * mm;
  physics_list->SetDefaultCutValue(production_cut);
  G4cout << "Production cut: " << production_cut / mm << " mm" << G4endl;
  run_manager->SetUserInitialization(physics_list);

  // Optical parameters (Cherenkov)
  auto optical_params = G4OpticalParameters::Instance();
  optical_params->SetCerenkovMaxPhotonsPerStep(100);
  optical_params->SetCerenkovStackPhotons(true);
  optical_params->SetCerenkovTrackSecondariesFirst(true);
  optical_params->SetCerenkovVerboseLevel(1);
  optical_params->SetBoundaryVerboseLevel(1);
  optical_params->SetAbsorptionVerboseLevel(1);

  run_manager->SetUserInitialization(new ActionInitialization());
  run_manager->Initialize();

  G4VisManager* vis_manager = new G4VisExecutive("Quiet");
  vis_manager->Initialize();

  G4UImanager* ui_manager = G4UImanager::GetUIpointer();

  if (!macro.empty()) {
    ui_manager->ApplyCommand("/control/execute " + macro);
  } else {
    ui_manager->ApplyCommand("/control/execute vis.mac");
    if (ui->IsGUI())
      ui_manager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }

  gAnaMan.CloseOutputFile();

  delete vis_manager;
  delete run_manager;

  return 0;
}
