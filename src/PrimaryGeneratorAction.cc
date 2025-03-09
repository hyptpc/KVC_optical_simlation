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

namespace
{
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::GeV;
  const auto particleTable = G4ParticleTable::GetParticleTable();
  auto& gAnaMan = AnaManager::GetInstance();
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun(1); // Initialize particle gun
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  GenerateKaon(anEvent);
}

void PrimaryGeneratorAction::GenerateKaon(G4Event* anEvent)
{
  static const G4String particle_name = "kaon-";
  // static const G4String particle_name = "pi-";
  static const auto particle = particleTable->FindParticle(particle_name);
  fParticleGun->SetParticleDefinition(particle);

  
  // -----------------------
  // 🎯 ビームエネルギーを正規分布で振る
  // -----------------------
  G4double E0 = 0.735 * GeV;
  G4double sigmaE = E0 * 0.02 / 2.355;
  G4double energy = G4RandGauss::shoot(E0, sigmaE);
  gAnaMan.SetBeamEnergy(energy);

  // 運動量を計算 (p = sqrt(E^2 - m^2))
  G4double mass = particle->GetPDGMass();
  G4double momentum = std::sqrt(energy * energy - mass * mass);

  // -----------------------
  // 📍 初期位置を正規分布で振る
  // -----------------------
  G4double x0 = 0.0 * mm, sigmaX = 50.0 * mm;
  G4double y0 = 0.0 * mm, sigmaY = 50.0 * mm;
  G4double z0 = -100.0 * mm;  // ビームの初期位置 (zは固定)

  // G4double x = G4RandGauss::shoot(x0, sigmaX);
  // G4double y = G4RandGauss::shoot(y0, sigmaY);
  G4double x = 0.0 * mm;
  G4double y = 0.0 * mm;
  G4double z = z0;

  G4ThreeVector position(x, y, z);
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  // -----------------------
  // 🎯 ビームの角度をランダムに振る
  // -----------------------
  // G4double theta_max = 0.1 * deg;  // 角度範囲（最大5度）
  // G4double theta = G4UniformRand() * theta_max;  // 0〜5度の一様乱数
  // G4double phi = G4UniformRand() * 360.0 * deg;  // 0〜360度の一様乱数
  G4double theta = 0.;
  G4double phi = 0.;
  
  G4double px = momentum * std::sin(theta) * std::cos(phi);
  G4double py = momentum * std::sin(theta) * std::sin(phi);
  G4double pz = momentum * std::cos(theta);

  G4ThreeVector direction(px, py, pz);
  gAnaMan.SetBeamMomentum(direction);
  direction = direction.unit();  // 正規化

  fParticleGun->SetParticleMomentum(momentum);
  fParticleGun->SetParticleMomentumDirection(direction);
  
  // -----------------------
  // 🛠 デバッグ用出力
  // -----------------------
  G4cout << "Particle: " << particle->GetParticleName()
	 << " | Energy: " << energy / GeV << " GeV"
	 << " | Momentum: " << momentum / GeV << " GeV/c"
	 << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm"
	 << " | Direction: (" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")"
	 << G4endl;

  // -----------------------
  // 🎯 粒子を生成
  // -----------------------
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

