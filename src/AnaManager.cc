// -*- C++ -*-

#include "AnaManager.hh"

#include <G4Event.hh>
#include <G4HCofThisEvent.hh>
#include <G4Run.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>

#include <TFile.h>
#include <TTree.h>

#include "MPPCHit.hh"

#define DEBUG 0

//_____________________________________________________________________________
AnaManager&
AnaManager::GetInstance()
{
  static AnaManager s_instance;
  return s_instance;
}

//_____________________________________________________________________________
AnaManager::AnaManager()
  : m_output_rootfile_path("test.root"),
    m_file(nullptr),
    m_tree(nullptr),
    m_evnum(0),
    m_event_id(0),
    m_nhit_mppc(0),
    m_cerenkov_all(0),
    m_cerenkov_quartz(0),
    m_beam_energy(0.),
    m_beam_mom_x(0.),
    m_beam_mom_y(0.),
    m_beam_mom_z(0.),
    m_beam_pos_x(0.),
    m_beam_pos_y(0.),
    m_beam_pos_z(0.),
    m_n_cherenkov_gen(0),
    m_n_delta_electrons(0),
    m_npe(0),
    m_n_trapped_air(0)
{
}

//_____________________________________________________________________________
AnaManager::~AnaManager()
{
}

//_____________________________________________________________________________
void
AnaManager::BeginOfRunAction(const G4Run* /* aRun */)
{
  // The output file and tree are created only once, at the first run.
  // Events of all runs in one session are stored in the same tree.
  if (m_file) return;

  m_file = new TFile(m_output_rootfile_path, "RECREATE");
  // Create the tree inside the output file so that baskets are written to disk
  m_tree = new TTree("tree", "GEANT4 optical simulation for KVC");

  m_tree->Branch("evnum", &m_evnum, "evnum/I");
  m_tree->Branch("event_id", &m_event_id, "event_id/I");
  m_tree->Branch("cerenkov_all", &m_cerenkov_all, "cerenkov_all/I");
  m_tree->Branch("cerenkov_quartz", &m_cerenkov_quartz, "cerenkov_quartz/I");

  // Beam info
  m_tree->Branch("beam_energy", &m_beam_energy, "beam_energy/D");
  m_tree->Branch("beam_mom_x", &m_beam_mom_x, "beam_mom_x/D");
  m_tree->Branch("beam_mom_y", &m_beam_mom_y, "beam_mom_y/D");
  m_tree->Branch("beam_mom_z", &m_beam_mom_z, "beam_mom_z/D");
  m_tree->Branch("beam_pos_x", &m_beam_pos_x, "beam_pos_x/D");
  m_tree->Branch("beam_pos_y", &m_beam_pos_y, "beam_pos_y/D");
  m_tree->Branch("beam_pos_z", &m_beam_pos_z, "beam_pos_z/D");
  m_tree->Branch("n_cherenkov_gen", &m_n_cherenkov_gen, "n_cherenkov_gen/I"); // Generated Cherenkov photons
  m_tree->Branch("n_delta_e", &m_n_delta_electrons, "n_delta_e/I");           // Generated delta electrons
  m_tree->Branch("npe", &m_npe, "npe/I");                                     // Detected photoelectrons

  // Trapping / monitoring info
  m_tree->Branch("nTrapped_Air", &m_n_trapped_air, "nTrapped_Air/I");

  // MPPC info
  m_tree->Branch("nhit_mppc", &m_nhit_mppc, "nhit_mppc/I");
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
void
AnaManager::BeginOfEventAction(const G4Event* /* anEvent */)
{
  m_n_trapped_air = 0;
  m_gen_wave_length.clear();
}

//_____________________________________________________________________________
void
AnaManager::EndOfEventAction(const G4Event* anEvent)
{
  m_event_id = anEvent->GetEventID();

  m_nhit_mppc = 0;
  m_npe = 0;

  // A missing hits collection is treated as zero hits so that every event is filled
  G4THitsCollection<MPPCHit>* mppc_hc = nullptr;
  G4HCofThisEvent* HCTE = anEvent->GetHCofThisEvent();
  const G4int mppc_hc_id = G4SDManager::GetSDMpointer()->GetCollectionID("MppcCollection");
  if (HCTE && mppc_hc_id >= 0) {
    mppc_hc = dynamic_cast<G4THitsCollection<MPPCHit>*>(HCTE->GetHC(mppc_hc_id));
    if (mppc_hc) {
      m_nhit_mppc = mppc_hc->entries();
    }
  }

  ResetContainer();
  for (G4int i = 0; i < m_nhit_mppc; ++i) {
    const MPPCHit* hit = (*mppc_hc)[i];

    const G4ThreeVector pos = hit->GetPosition();
    m_pos_x.push_back(pos.x());
    m_pos_y.push_back(pos.y());
    m_pos_z.push_back(pos.z());
    m_time.push_back(hit->GetTime());
    m_energy.push_back(hit->GetEnergy());
    m_wave_length.push_back(hit->GetWaveLength());
    m_particle_id.push_back(hit->GetParticleID());
    m_seg.push_back(hit->GetCopyNumber());

    const G4int detect_flag = hit->GetDetectFlag();
    m_detect_flag.push_back(detect_flag);
    if (detect_flag == 1) ++m_npe;
  }

  m_tree->Fill();
  ++m_evnum;
#if DEBUG
  G4cout << m_evnum << ", " << m_nhit_mppc << G4endl;
#endif
}

//_____________________________________________________________________________
void
AnaManager::EndOfRunAction(const G4Run* /* aRun */)
{
  // Write the tree after each run so that the file is valid even if the session ends abnormally
  if (m_file && m_file->IsOpen()) {
    m_file->cd();
    m_tree->Write("", TObject::kOverwrite);
  }
}

//_____________________________________________________________________________
void
AnaManager::CloseOutputFile()
{
  if (!m_file) return;
  if (m_file->IsOpen()) {
    m_file->cd();
    m_tree->Write("", TObject::kOverwrite);
    m_file->Close();
  }
  delete m_file; // also deletes the tree owned by the file
  m_file = nullptr;
  m_tree = nullptr;
}

//_____________________________________________________________________________
void
AnaManager::ResetContainer()
{
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

//_____________________________________________________________________________
void
AnaManager::SetNumOfCerenkovAll(G4int cerenkov_all)
{
  m_cerenkov_all = cerenkov_all;
}

//_____________________________________________________________________________
void
AnaManager::SetNumOfCerenkovQuartz(G4int cerenkov_quartz)
{
  m_cerenkov_quartz = cerenkov_quartz;
}

//_____________________________________________________________________________
void
AnaManager::SetBeamEnergy(G4double beam_energy)
{
  m_beam_energy = beam_energy;
}

//_____________________________________________________________________________
void
AnaManager::SetBeamMomentum(const G4ThreeVector& beam_momentum)
{
  m_beam_mom_x = beam_momentum.x();
  m_beam_mom_y = beam_momentum.y();
  m_beam_mom_z = beam_momentum.z();
}

//_____________________________________________________________________________
void
AnaManager::SetBeamPosition(const G4ThreeVector& beam_position)
{
  m_beam_pos_x = beam_position.x();
  m_beam_pos_y = beam_position.y();
  m_beam_pos_z = beam_position.z();
}

//_____________________________________________________________________________
void
AnaManager::SetOutputRootfilePath(const G4String& output_rootfile_path)
{
  m_output_rootfile_path = output_rootfile_path;
}

//_____________________________________________________________________________
G4String
AnaManager::GetOutputRootfilePath() const
{
  return m_output_rootfile_path;
}
