#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"

// Forward declarations for ROOT classes
class TFile;
class TTree;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
  G4ParticleGun* fParticleGun; // Particle gun
  void GenerateBeam(G4Event* anEvent);
  void GeneratePhoton(G4Event* anEvent);
  void GenerateRootBeam(G4Event* anEvent);

  // ROOT beam members
  TFile* fRootFile;
  TTree* fTree;
  G4int  fMaxEntries;

  // Branch variables
  G4double fPx, fPy, fPz;
  G4double fVx, fVy, fVz;
};

#endif
