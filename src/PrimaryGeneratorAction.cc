#include "PrimaryGeneratorAction.hh"
#include "AnaManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

#include <TFile.h>
#include <TTree.h>

#include "ConfManager.hh"

#define DEBUG 0

namespace
{
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::GeV;
  const auto particleTable = G4ParticleTable::GetParticleTable();
  auto& gAnaMan  = AnaManager::GetInstance();
  auto& gConfMan = ConfManager::GetInstance();
}

//_____________________________________________________________________________
PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fRootFile(nullptr), fTree(nullptr), fMaxEntries(0)
{
  fParticleGun = new G4ParticleGun(1);

  // Initialize ROOT beam if file is provided
  G4String input_file = gConfMan.Get("input_beam_file");
  if(!input_file.empty() && input_file != "none") {
    fRootFile = new TFile(input_file, "READ");
    if(fRootFile && fRootFile->IsOpen()) {
      fTree = (TTree*)fRootFile->Get("tree"); // Expecting tree named "tree"
      if(fTree) {
        fMaxEntries = fTree->GetEntries();
        fTree->SetBranchAddress("px", &fPx);
        fTree->SetBranchAddress("py", &fPy);
        fTree->SetBranchAddress("pz", &fPz);
        fTree->SetBranchAddress("vx", &fVx);
        fTree->SetBranchAddress("vy", &fVy);
        fTree->SetBranchAddress("vz", &fVz);
        G4cout << "PrimaryGeneratorAction: ROOT beam mode enabled (Random Sampling). File: " << input_file 
               << " (Pool Size: " << fMaxEntries << ")" << G4endl;
      } else {
        G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction", "TreeNotFound", FatalException, "TTree 'tree' not found in ROOT beam file.");
      }
    } else {
      G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction", "FileNotFound", FatalException, "Failed to open ROOT beam file.");
    }
  }
}

//_____________________________________________________________________________
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  if(fRootFile) {
    fRootFile->Close();
    delete fRootFile;
  }
}

//_____________________________________________________________________________
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fTree) {
    GenerateRootBeam(anEvent);
  } else {
    GenerateBeam(anEvent);
  }
}

//_____________________________________________________________________________
void PrimaryGeneratorAction::GenerateBeam(G4Event *anEvent)
{
  static const G4String particle_name = gConfMan.Get("particle");
  static const auto particle = particleTable->FindParticle(particle_name);
  fParticleGun->SetParticleDefinition(particle);

  // -----------------------
  // Momentum
  // -----------------------
  G4double p0 = gConfMan.GetDouble("momentum") * GeV;
  G4double sigma_p = p0 * 0.02 / 2.355;
  // G4double momentum = G4RandGauss::shoot(p0, sigma_p);
  G4double momentum = p0;

  G4double mass = particle->GetPDGMass();
  G4double energy = std::sqrt(mass * mass + momentum * momentum);
  G4double kineticE = energy - mass;

  gAnaMan.SetBeamEnergy(kineticE);
  fParticleGun->SetParticleEnergy(kineticE); // Set kinetic energy

  // -----------------------
  // Momentum direction
  // -----------------------
  // G4double theta_max = 0.1 * deg;
  // G4double theta = G4UniformRand() * theta_max;
  // G4double phi = G4UniformRand() * 360.0 * deg;
  G4double theta = 0.;
  G4double phi = 0.;

  G4double px = momentum * std::sin(theta) * std::cos(phi);
  G4double py = momentum * std::sin(theta) * std::sin(phi);
  G4double pz = momentum * std::cos(theta);

  G4ThreeVector direction(px, py, pz);
  gAnaMan.SetBeamMomentum(direction);
  direction = direction.unit(); // normalize

  // fParticleGun->SetParticleMomentum(momentum);
  fParticleGun->SetParticleMomentumDirection(direction);

  // -----------------------
  // Position
  // -----------------------
  // G4double x0 = 0.0 * mm, sigmaX = 1.0 * mm;
  // G4double y0 = 0.0 * mm, sigmaY = 1.0 * mm;
  G4double z0 = -100.0 * mm;

  // G4double x = G4RandGauss::shoot(x0, sigmaX);
  // G4double y = G4RandGauss::shoot(y0, sigmaY);
  G4double x = 0.0 * mm;
  G4double y = gConfMan.GetDouble("beam_y_offset") * mm;
  G4double z = z0;

  G4ThreeVector position(x, y, z);
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  // -----------------------
  // Debug
  // -----------------------
#if DEBUG
  G4cout << "Particle: " << particle->GetParticleName() << G4endl
  	 << " | Energy: " << energy / GeV << " GeV" << G4endl
  	 << " | Momentum: " << momentum / GeV << " GeV/c" << G4endl
  	 << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm"  << G4endl
  	 << " | Direction: (" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")"
  	 << G4endl;
#endif

  // Gun
  // -----------------------
  fParticleGun->GeneratePrimaryVertex(anEvent);
}


//_____________________________________________________________________________
void PrimaryGeneratorAction::GeneratePhoton(G4Event* anEvent)
{
  static const G4String particle_name = "opticalphoton";
  static const auto particle = particleTable->FindParticle(particle_name);
  fParticleGun->SetParticleDefinition(particle);

  
  // -----------------------
  // Energy
  // -----------------------
  G4double wl_min = 320. * CLHEP::nm;
  G4double wl_max = 900. * CLHEP::nm;
  
  G4double wavelength = G4UniformRand() * (wl_max - wl_min) + wl_min;
  wavelength = 400.0 * CLHEP::nm;
  G4double energy = (CLHEP::h_Planck * CLHEP::c_light  / wavelength);
  gAnaMan.SetBeamEnergy(energy);  
  fParticleGun->SetParticleEnergy(energy);

    
  // -----------------------
  // direction
  // -----------------------
  // G4double beta_min = 0.794;
  // G4double beta_max = 0.847;
  G4double beta_min = 0.95;
  G4double beta_max = 1.0;
  
  G4double beta = G4UniformRand() * (beta_max - beta_min) + beta_min;
  beta = 0.83;
  G4double theta = std::acos(1./(1.46*beta));
  G4double phi = G4UniformRand() * 360.0 * deg;
  phi = 0.0;
  G4double px = std::sin(theta) * std::cos(phi);
  G4double py = std::sin(theta) * std::sin(phi);
  G4double pz = std::cos(theta);

  G4ThreeVector direction(px, py, pz);
  gAnaMan.SetBeamMomentum(direction);
  direction = direction.unit();  // normalize

  // fParticleGun->SetParticleMomentum(momentum);
  fParticleGun->SetParticleMomentumDirection(direction);

  // -----------------------
  // Position
  // -----------------------
  G4double x0 = 0.0 * mm, sigmaX = 1.0 * mm;
  G4double y0 = 0.0 * mm, sigmaY = 1.0 * mm;
  G4double z0 = (G4UniformRand() * 18.0 - 11.0) * mm;
  z0 = 0.0;
  // G4double x = G4RandGauss::shoot(x0, sigmaX);
  // G4double y = G4RandGauss::shoot(y0, sigmaY);
  G4double x = 0.0 * mm;
  G4double y = 0.0 * mm;
  G4double z = z0;

  G4ThreeVector position(x, y, z);
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  // -----------------------
  // Debug
  // -----------------------
#if DEBUG
  G4double momentum = energy;
  G4cout << "Particle: " << particle->GetParticleName() << G4endl
	 << " | Energy: " << energy / CLHEP::eV << " eV" << G4endl
	 << " | Momentum: " << momentum / GeV << " GeV/c" << G4endl
	 << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm"  << G4endl
	 << " | Direction: (" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")"
	 << G4endl;
#endif

  // Gun
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//_____________________________________________________________________________
void PrimaryGeneratorAction::GenerateRootBeam(G4Event* anEvent)
{
  // Random sampling (Bootstrap) from the pool of realistic particles
  G4int entry = G4RandFlat::shootInt(fMaxEntries);
  fTree->GetEntry(entry);

  static const G4String particle_name = gConfMan.Get("particle");
  static const auto particle = particleTable->FindParticle(particle_name);
  fParticleGun->SetParticleDefinition(particle);

  G4ThreeVector direction(fPx, fPy, fPz);
  G4double momentum = direction.mag() * MeV; // Fixed: Input is MeV/c
  direction = direction.unit();

  G4double mass = particle->GetPDGMass();
  G4double energy = std::sqrt(mass * mass + momentum * momentum);
  G4double kineticE = energy - mass;

  gAnaMan.SetBeamEnergy(kineticE);
  gAnaMan.SetBeamMomentum(direction * momentum);
  fParticleGun->SetParticleEnergy(kineticE);
  fParticleGun->SetParticleMomentumDirection(direction);

  G4double thickness = gConfMan.GetDouble("quartz_thickness") * mm;
  G4double z_surf = -thickness / 2.0;

  // ROOT file Z is ~ -10 mm (relative to surface). 
  // We align this to Geant4 surface position.
  G4double y_offset = gConfMan.GetDouble("beam_y_offset") * mm;
  G4ThreeVector position(fVx * mm, fVy * mm + y_offset, z_surf + fVz * mm);
  
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
