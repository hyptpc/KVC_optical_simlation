// -*- C++ -*-

#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include <G4Types.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

class G4Event;
class G4ParticleGun;
class TFile;
class TTree;

//_____________________________________________________________________________
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

public:
  void GeneratePrimaries(G4Event* anEvent) override;

private:
  void GenerateBeam(G4Event* anEvent);     // Particle gun (conf: particle, momentum)
  void GeneratePhoton(G4Event* anEvent);   // Single optical photon (for tests)
  void GenerateRootBeam(G4Event* anEvent); // Sampled from the beam file

private:
  G4ParticleGun* m_particle_gun;

  // ROOT beam file (conf: input_beam_file)
  TFile* m_beam_file;
  TTree* m_beam_tree;
  G4int  m_n_beam_entries;

  // Branch variables of the beam tree
  G4double m_px, m_py, m_pz; // [MeV/c]
  G4double m_vx, m_vy, m_vz; // [mm]
};

#endif
