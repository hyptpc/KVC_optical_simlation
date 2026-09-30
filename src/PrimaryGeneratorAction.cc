// -*- C++ -*-

#include "PrimaryGeneratorAction.hh"

#include <cmath>

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TTree.h>

#include "AnaManager.hh"
#include "ConfManager.hh"

#define DEBUG 0

namespace
{
  using CLHEP::deg;
  using CLHEP::GeV;
  using CLHEP::mm;
  const auto particle_table = G4ParticleTable::GetParticleTable();
  auto& gAnaMan  = AnaManager::GetInstance();
  auto& gConfMan = ConfManager::GetInstance();
}

//_____________________________________________________________________________
PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    m_particle_gun(new G4ParticleGun(1)),
    m_beam_file(nullptr),
    m_beam_tree(nullptr),
    m_n_beam_entries(0),
    m_px(0.), m_py(0.), m_pz(0.),
    m_vx(0.), m_vy(0.), m_vz(0.)
{
  // Initialize ROOT beam if file is provided
  const G4String input_file = gConfMan.Get("input_beam_file");
  if (input_file.empty() || input_file == "none") return;

  // Relative path is resolved against the conf file directory
  const G4String beam_file_path = gConfMan.GetPath("input_beam_file");
  m_beam_file = new TFile(beam_file_path, "READ");
  if (!m_beam_file->IsOpen()) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction", "FileNotFound",
                FatalException, "Failed to open ROOT beam file.");
    return;
  }

  m_beam_tree = dynamic_cast<TTree*>(m_beam_file->Get("tree")); // Expecting tree named "tree"
  if (!m_beam_tree) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction", "TreeNotFound",
                FatalException, "TTree 'tree' not found in ROOT beam file.");
    return;
  }

  m_n_beam_entries = m_beam_tree->GetEntries();
  m_beam_tree->SetBranchAddress("px", &m_px);
  m_beam_tree->SetBranchAddress("py", &m_py);
  m_beam_tree->SetBranchAddress("pz", &m_pz);
  m_beam_tree->SetBranchAddress("vx", &m_vx);
  m_beam_tree->SetBranchAddress("vy", &m_vy);
  m_beam_tree->SetBranchAddress("vz", &m_vz);
  G4cout << "PrimaryGeneratorAction: ROOT beam mode enabled (Random Sampling). File: "
         << beam_file_path << " (Pool Size: " << m_n_beam_entries << ")" << G4endl;
}

//_____________________________________________________________________________
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete m_particle_gun;
  if (m_beam_file) {
    m_beam_file->Close();
    delete m_beam_file;
  }
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (m_beam_tree) {
    GenerateRootBeam(anEvent);
  } else {
    GenerateBeam(anEvent);
  }
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateBeam(G4Event* anEvent)
{
  static const G4String particle_name = gConfMan.Get("particle");
  static const auto particle = particle_table->FindParticle(particle_name);
  m_particle_gun->SetParticleDefinition(particle);

  // -----------------------
  // Momentum
  // -----------------------
  const G4double p0 = gConfMan.GetDouble("momentum") * GeV;
  // const G4double sigma_p = p0 * 0.02 / 2.355;
  // const G4double momentum = G4RandGauss::shoot(p0, sigma_p);
  const G4double momentum = p0;

  const G4double mass = particle->GetPDGMass();
  const G4double energy = std::sqrt(mass * mass + momentum * momentum);
  const G4double kinetic_energy = energy - mass;

  gAnaMan.SetBeamEnergy(kinetic_energy);
  m_particle_gun->SetParticleEnergy(kinetic_energy);

  // -----------------------
  // Momentum direction
  // -----------------------
  // const G4double theta_max = 0.1 * deg;
  // const G4double theta = G4UniformRand() * theta_max;
  // const G4double phi = G4UniformRand() * 360.0 * deg;
  const G4double theta = 0.;
  const G4double phi = 0.;

  const G4double px = momentum * std::sin(theta) * std::cos(phi);
  const G4double py = momentum * std::sin(theta) * std::sin(phi);
  const G4double pz = momentum * std::cos(theta);

  G4ThreeVector direction(px, py, pz);
  gAnaMan.SetBeamMomentum(direction);
  direction = direction.unit();
  m_particle_gun->SetParticleMomentumDirection(direction);

  // -----------------------
  // Position
  // -----------------------
  // const G4double x0 = 0.0 * mm, sigma_x = 1.0 * mm;
  // const G4double y0 = 0.0 * mm, sigma_y = 1.0 * mm;
  // const G4double x = G4RandGauss::shoot(x0, sigma_x);
  // const G4double y = G4RandGauss::shoot(y0, sigma_y);
  const G4double x = 0.0 * mm;
  const G4double y = gConfMan.GetDouble("beam_y_offset") * mm;
  const G4double z = -100.0 * mm;

  const G4ThreeVector position(x, y, z);
  m_particle_gun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

#if DEBUG
  G4cout << "Particle: " << particle->GetParticleName() << G4endl
         << " | Energy: " << energy / GeV << " GeV" << G4endl
         << " | Momentum: " << momentum / GeV << " GeV/c" << G4endl
         << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm" << G4endl
         << " | Direction: (" << direction.x() << ", " << direction.y() << ", "
         << direction.z() << ")" << G4endl;
#endif

  m_particle_gun->GeneratePrimaryVertex(anEvent);
}

//_____________________________________________________________________________
// Shoot a single optical photon (for tests; currently not used)
void
PrimaryGeneratorAction::GeneratePhoton(G4Event* anEvent)
{
  static const G4String particle_name = "opticalphoton";
  static const auto particle = particle_table->FindParticle(particle_name);
  m_particle_gun->SetParticleDefinition(particle);

  // -----------------------
  // Energy
  // -----------------------
  // const G4double wl_min = 320. * CLHEP::nm;
  // const G4double wl_max = 900. * CLHEP::nm;
  // const G4double wave_length = G4UniformRand() * (wl_max - wl_min) + wl_min;
  const G4double wave_length = 400.0 * CLHEP::nm;
  const G4double energy = (CLHEP::h_Planck * CLHEP::c_light / wave_length);
  gAnaMan.SetBeamEnergy(energy);
  m_particle_gun->SetParticleEnergy(energy);

  // -----------------------
  // Direction (Cherenkov angle in quartz for a given beta)
  // -----------------------
  // const G4double beta_min = 0.95;
  // const G4double beta_max = 1.0;
  // const G4double beta = G4UniformRand() * (beta_max - beta_min) + beta_min;
  const G4double beta = 0.83;
  const G4double theta = std::acos(1. / (1.46 * beta));
  // const G4double phi = G4UniformRand() * 360.0 * deg;
  const G4double phi = 0.0;
  const G4double px = std::sin(theta) * std::cos(phi);
  const G4double py = std::sin(theta) * std::sin(phi);
  const G4double pz = std::cos(theta);

  G4ThreeVector direction(px, py, pz);
  gAnaMan.SetBeamMomentum(direction);
  direction = direction.unit();
  m_particle_gun->SetParticleMomentumDirection(direction);

  // -----------------------
  // Position
  // -----------------------
  // const G4double z = (G4UniformRand() * 18.0 - 11.0) * mm;
  const G4double x = 0.0 * mm;
  const G4double y = 0.0 * mm;
  const G4double z = 0.0 * mm;

  const G4ThreeVector position(x, y, z);
  m_particle_gun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

#if DEBUG
  G4cout << "Particle: " << particle->GetParticleName() << G4endl
         << " | Energy: " << energy / CLHEP::eV << " eV" << G4endl
         << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm" << G4endl
         << " | Direction: (" << direction.x() << ", " << direction.y() << ", "
         << direction.z() << ")" << G4endl;
#endif

  m_particle_gun->GeneratePrimaryVertex(anEvent);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateRootBeam(G4Event* anEvent)
{
  // Random sampling (bootstrap) from the pool of measured beam particles
  const G4int entry = G4RandFlat::shootInt(m_n_beam_entries);
  m_beam_tree->GetEntry(entry);

  static const G4String particle_name = gConfMan.Get("particle");
  static const auto particle = particle_table->FindParticle(particle_name);
  m_particle_gun->SetParticleDefinition(particle);

  G4ThreeVector direction(m_px, m_py, m_pz);
  const G4double momentum = direction.mag() * MeV; // Input is in MeV/c
  direction = direction.unit();

  const G4double mass = particle->GetPDGMass();
  const G4double energy = std::sqrt(mass * mass + momentum * momentum);
  const G4double kinetic_energy = energy - mass;

  gAnaMan.SetBeamEnergy(kinetic_energy);
  gAnaMan.SetBeamMomentum(direction * momentum);
  m_particle_gun->SetParticleEnergy(kinetic_energy);
  m_particle_gun->SetParticleMomentumDirection(direction);

  // The z position in the beam file is relative to the upstream surface of the quartz
  // (~ -10 mm). Align it to the surface position in Geant4.
  const G4double thickness = gConfMan.GetDouble("quartz_thickness") * mm;
  const G4double z_surface = -thickness / 2.0;
  const G4double y_offset = gConfMan.GetDouble("beam_y_offset") * mm;
  const G4ThreeVector position(m_vx * mm, m_vy * mm + y_offset, z_surface + m_vz * mm);

  m_particle_gun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  m_particle_gun->GeneratePrimaryVertex(anEvent);
}
