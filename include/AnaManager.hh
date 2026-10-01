// -*- C++ -*-

#ifndef ANA_MANAGER_HH
#define ANA_MANAGER_HH

#include <vector>

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>

class G4Event;
class G4Run;
class TFile;
class TTree;

//_____________________________________________________________________________
class AnaManager
{
public:
  static AnaManager& GetInstance();
  ~AnaManager();

private:
  AnaManager();
  AnaManager(const AnaManager&);
  AnaManager& operator=(const AnaManager&);

private:
  G4String m_output_rootfile_path;
  TFile*   m_file;
  TTree*   m_tree;

  // Event info
  G4int m_evnum;
  G4int m_event_id;
  G4int m_nhit_mppc;
  G4int m_cerenkov_all;
  G4int m_cerenkov_quartz;

  // Beam info
  G4double m_beam_energy;
  G4double m_beam_mom_x;
  G4double m_beam_mom_y;
  G4double m_beam_mom_z;
  G4double m_beam_pos_x;
  G4double m_beam_pos_y;
  G4double m_beam_pos_z;

  G4int m_n_cherenkov_gen;   // Number of generated Cherenkov photons
  G4int m_n_delta_electrons; // Number of generated delta electrons
  G4int m_npe;               // Number of detected photoelectrons
  G4int m_n_trapped_air;     // Number of photons from the quartz lost outside the quartz (not at the MPPCs)

  // Per-photon info
  std::vector<G4double> m_gen_wave_length; // Wavelengths of generated Cherenkov photons
  std::vector<G4double> m_pos_x;
  std::vector<G4double> m_pos_y;
  std::vector<G4double> m_pos_z;
  std::vector<G4double> m_time;
  std::vector<G4double> m_energy;
  std::vector<G4double> m_wave_length;
  std::vector<G4int>    m_particle_id;
  std::vector<G4int>    m_seg;
  std::vector<G4int>    m_detect_flag;

public:
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  void BeginOfEventAction(const G4Event* anEvent);
  void EndOfEventAction(const G4Event* anEvent);
  void CloseOutputFile(); // Write the tree and close the output file (call once at the end)

  void     ResetContainer();
  void     SetNumOfCerenkovAll(G4int cerenkov_all);
  void     SetNumOfCerenkovQuartz(G4int cerenkov_quartz);
  void     SetBeamEnergy(G4double beam_energy);
  void     SetBeamMomentum(const G4ThreeVector& beam_momentum);
  void     SetBeamPosition(const G4ThreeVector& beam_position);
  void     SetOutputRootfilePath(const G4String& output_rootfile_path);
  G4String GetOutputRootfilePath() const;
  void     SetCherenkovGen(G4int n) { m_n_cherenkov_gen = n; }
  void     SetNumDeltaElectrons(G4int n) { m_n_delta_electrons = n; }
  void     IncrementTrappedAir() { ++m_n_trapped_air; }
  void     AddGenWavelength(G4double wave_length) { m_gen_wave_length.push_back(wave_length); }
};

#endif
