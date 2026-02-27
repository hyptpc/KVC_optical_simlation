#include "AnaManager.hh"
#include "ConfManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"

#include "MPPCHit.hh"

#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include <string>
#include <sstream>
#include <vector>

#include "G4ThreeVector.hh"

extern int gCerenkovCounter;
extern double decay_check;

#define DEBUG 0

//_____________________________________________________________________________
AnaManager& AnaManager::GetInstance()
{
  static AnaManager instance;
  return instance;
}

AnaManager::AnaManager()
  : m_file(),
    m_output_rootfile_path("test.root"),
    m_tree(new TTree("tree", "GEANT4 optical simulation for KVC")),
    m_evnum(0),
    m_event_id(0),
    m_nhit_mppc(0),
    m_cerenkov_all(0),
    m_cerenkov_quartz(0),
    n_cherenkov_gen(0), // Number of generated Cherenkov photons
    m_beam_energy(0.),
    m_beam_mom_x(0.),
    m_beam_mom_y(0.),
    m_beam_mom_z(0.),
    m_beam_pos_x(0.),
    m_beam_pos_y(0.),
    m_beam_pos_z(0.),
    m_npe(0),          // Number of detected photoelectrons
    m_nTrapped_Air(0)
{
}

AnaManager::~AnaManager()
{
}


//_____________________________________________________________________________
void AnaManager::BeginOfRunAction(const G4Run*)
{
  m_file = new TFile(m_output_rootfile_path, "RECREATE");
  m_tree->Reset();

  m_tree->Branch("evnum", &m_evnum, "evnum/I");
  m_tree->Branch("event_id", &m_event_id, "event_id/I");
  m_tree->Branch("cerenkov_all", &m_cerenkov_all, "cerenkov_all/I");
  m_tree->Branch("cerenkov_quartz", &m_cerenkov_quartz, "cerenkov_quartz/I");

  // beam info
  m_tree->Branch("beam_energy", &m_beam_energy, "beam_energy/D");
  m_tree->Branch("beam_mom_x", &m_beam_mom_x, "beam_mom_x/D");
  m_tree->Branch("beam_mom_y", &m_beam_mom_y, "beam_mom_y/D");
  m_tree->Branch("beam_mom_z", &m_beam_mom_z, "beam_mom_z/D");
  m_tree->Branch("beam_pos_x", &m_beam_pos_x, "beam_pos_x/D");
  m_tree->Branch("beam_pos_y", &m_beam_pos_y, "beam_pos_y/D");
  m_tree->Branch("beam_pos_z", &m_beam_pos_z, "beam_pos_z/D");
  m_tree->Branch("n_cherenkov_gen", &n_cherenkov_gen, "n_cherenkov_gen/I"); // Number of generated Cherenkov photons
  m_tree->Branch("npe", &m_npe, "npe/I");           // Number of detected photoelectrons
  
  // Trapping/Monitoring info
  m_tree->Branch("nTrapped_Air",    &m_nTrapped_Air,    "nTrapped_Air/I");
  
  
  // MPPC info
  m_tree->Branch("nhit_mppc",&m_nhit_mppc,"nhit_mppc/I");
  m_tree->Branch("pos_x", &m_pos_x);
  m_tree->Branch("pos_y", &m_pos_y);
  m_tree->Branch("pos_z", &m_pos_z);
  m_tree->Branch("time", &m_time);
  m_tree->Branch("energy", &m_energy);
  m_tree->Branch("wave_length", &m_wave_length);
  m_tree->Branch("particle_id", &m_particle_id);
  m_tree->Branch("seg", &m_seg);
  m_tree->Branch("detect_flag", &m_detect_flag);
  m_tree->Branch("gen_wave_length", &m_gen_wave_length);
}

//_____________________________________________________________________________
void AnaManager::BeginOfEventAction(const G4Event* anEvent)
{
  m_nTrapped_Air = 0;
}

//_____________________________________________________________________________
void AnaManager::EndOfEventAction(const G4Event* anEvent)
{
  G4HCofThisEvent* HCTE = anEvent->GetHCofThisEvent();
  if(!HCTE) return;
  m_event_id = anEvent->GetEventID();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  m_nhit_mppc = 0;  
  m_npe = 0; // initialization
  G4THitsCollection<MPPCHit>* MPPCHC;
  G4int ColIdMPPC = SDMan->GetCollectionID("MppcCollection");
  if (ColIdMPPC >= 0) {
    MPPCHC = dynamic_cast<G4THitsCollection<MPPCHit>*>(HCTE->GetHC(ColIdMPPC));
    if (MPPCHC) {
      m_nhit_mppc = MPPCHC->entries();
    }
  }

  ResetContainer();
  for (int i=0; i<m_nhit_mppc; i++) {
    MPPCHit* aHit = (*MPPCHC)[i];

    G4ThreeVector pos = aHit->GetPosition();
    // m_pos.push_back(TVector3(pos.x(), pos.y(), pos.z()));
    m_pos_x.push_back(pos.x());
    m_pos_y.push_back(pos.y());
    m_pos_z.push_back(pos.z());

    G4double time = aHit->GetTime();
    m_time.push_back(time);

    G4double energy = aHit->GetEnergy();
    m_energy.push_back(energy);

    G4double wave_length = aHit->GetWaveLength();
    m_wave_length.push_back(wave_length);
    
    G4int particle_id = aHit->GetParticleID();
    m_particle_id.push_back(particle_id);

    G4int seg = aHit->GetCopyNumber();
    m_seg.push_back(seg);

    G4int detect_flag = aHit->GetDetectFlag();
    m_detect_flag.push_back(detect_flag);
    if(detect_flag == 1) m_npe++; // count
  }
  
  m_tree->Fill();
  m_evnum++;
#if DEBUG
  G4cout << m_evnum << ", " << m_nhit_mppc << G4endl;
#endif
}

//_____________________________________________________________________________
void AnaManager::EndOfRunAction(const G4Run* aRun) {
  if (m_file && m_file->IsOpen()) {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
  }
}

//_____________________________________________________________________________
void AnaManager::ResetContainer()
{
  m_gen_wave_length.clear();
  m_pos_x.clear();
  m_pos_y.clear();
  m_pos_z.clear();
  m_time.clear();
  m_energy.clear();
  m_wave_length.clear();
  m_particle_id.clear();
  m_seg.clear();
  m_detect_flag.clear();
}

void AnaManager::SetNumOfCerenkovAll(G4int cerenkov_all)
{
  m_cerenkov_all = cerenkov_all;
}

void AnaManager::SetNumOfCerenkovQuartz(G4int cerenkov_quartz)
{
  m_cerenkov_quartz = cerenkov_quartz;
}

void AnaManager::SetBeamEnergy(G4double beam_energy)
{
  m_beam_energy = beam_energy;
}

void AnaManager::SetBeamMomentum(G4ThreeVector beam_momentum)
{
  m_beam_mom_x = beam_momentum.x();
  m_beam_mom_y = beam_momentum.y();
  m_beam_mom_z = beam_momentum.z();
}

void AnaManager::SetBeamPosition(G4ThreeVector beam_position)
{
  m_beam_pos_x = beam_position.x();
  m_beam_pos_y = beam_position.y();
  m_beam_pos_z = beam_position.z();
}

void AnaManager::SetOutputRootfilePath(G4String output_rootfile_path)
{
  m_output_rootfile_path = output_rootfile_path;
}
G4String AnaManager::GetOutputRootfilePath()
{
  return m_output_rootfile_path;
}
