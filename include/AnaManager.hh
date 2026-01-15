#ifndef ANA_MANAGER_HH
#define ANA_MANAGER_HH

#include <vector>

#include <G4ThreeVector.hh>
#include <G4Run.hh>
#include <G4Event.hh>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

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
  TFile* m_file;
  TTree* m_tree;
  G4int m_evnum;
  G4int m_nhit_mppc;
  G4int m_cerenkov_all;
  G4int m_cerenkov_quartz;
  G4double m_beam_energy;
  G4double m_beam_mom_x;
  G4double m_beam_mom_y;
  G4double m_beam_mom_z;
  G4double m_beam_pos_x;
  G4double m_beam_pos_y;
  G4double m_beam_pos_z;  
  G4int n_cer_gen; //チェレンコフ光発生数

  
  // std::vector<TVector3> m_pos;
  std::vector<G4double> m_pos_x;
  std::vector<G4double> m_pos_y;
  std::vector<G4double> m_pos_z;
  std::vector<G4double> m_time;
  std::vector<G4double> m_energy;
  std::vector<G4double> m_wave_length;
  std::vector<G4int> m_particle_id;
  std::vector<G4int> m_seg;
  std::vector<G4int> m_detect_flag;
    
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

  void ResetContainer();
  void SetNumOfCerenkovAll(G4int cerenkov_all);
  void SetNumOfCerenkovQuartz(G4int cerenkov_quartz);
  void SetBeamEnergy(G4double beam_energy);
  void SetBeamMomentum(G4ThreeVector beam_momentum);
  void SetBeamPosition(G4ThreeVector beam_position);
  void SetOutputRootfilePath(G4String output_rootfile_path);
  G4String GetOutputRootfilePath();
  //void FillTree(); //追加
  void SetCerGen(Int_t val) { n_cer_gen = val; } //追加

};

#endif
