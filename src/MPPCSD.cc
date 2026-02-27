#include "MPPCSD.hh"
#include "ConfManager.hh"
#include "KVC_OpticalProperties.hh"

#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "Randomize.hh"

#include "TGraph.h"
#include "TSpline.h"

//_____________________________________________________________________________
MPPCSD::MPPCSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_qe_spline(nullptr),
    m_range_min(1. * CLHEP::eV),
    m_range_max(7. * CLHEP::eV),
    m_qe_scale(1.0)
{
    collectionName.insert("MppcCollection");

    InitializeQESpline();

    m_qe_scale = ConfManager::GetInstance().GetDouble("qe_scale");
    if (m_qe_scale <= 0.0) m_qe_scale = 1.0;
}

//_____________________________________________________________________________
MPPCSD::~MPPCSD() {
  delete m_qe_spline;
}

//_____________________________________________________________________________
void MPPCSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<MPPCHit>(SensitiveDetectorName,
						     collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool MPPCSD::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
  const auto postStepPoint = aStep->GetPostStepPoint();  // step ends inside MPPC: use post for hit volume
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4int particleID = Definition->GetPDGEncoding();
  if (Definition != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

  G4ThreeVector worldPos = postStepPoint->GetPosition();
  G4ThreeVector pos      = postStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  G4double hitTime = postStepPoint->GetGlobalTime();
  G4double energy = aTrack->GetTotalEnergy();
  G4double waveLength = (CLHEP::h_Planck * CLHEP::c_light / energy) / CLHEP::nm;
  G4int copyNumber = postStepPoint->GetTouchableHandle()->GetCopyNumber();
  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  G4int detectFlag = 0;

  // -- kill track -----
  // NOTE: Optical photons entering MPPC are absorbed here regardless of QE result.
  // Photon detection is determined by detectFlag based on QE.
  aTrack->SetTrackStatus(fStopAndKill);

#ifdef USE_SURFACE_PDE
  // In Method A, Geant4 boundary process handled efficiency. 
  // Any photon reaching here is a survivor.
  detectFlag = 1;
#else
  // -- QE check (Method B) -----
  G4double eval_energy = energy;
  if      (eval_energy < m_range_min) eval_energy = m_range_min;
  else if (eval_energy > m_range_max) eval_energy = m_range_max;

  G4double qe_value = m_qe_spline->Eval(eval_energy) * m_qe_scale;
  if (qe_value > 1.0) qe_value = 1.0;

  G4double random_value = G4UniformRand();
  if (random_value <= qe_value) {
    detectFlag = 1;
  }
#endif

  // -- record -----
  if (detectFlag == 1) {
    MPPCHit* aHit = new MPPCHit();
    aHit->SetPosition(pos);
    aHit->SetWorldPosition(worldPos);
    aHit->SetEnergy(energy);
    aHit->SetWaveLength(waveLength);
    aHit->SetTime(hitTime);
    aHit->SetParticleID(particleID);
    aHit->SetCopyNumber(copyNumber);
    aHit->SetEventID(eventID);
    aHit->SetDetectFlag(detectFlag);

    m_hits_collection->insert(aHit);
  }

  
  return true;
}

//_____________________________________________________________________________
void MPPCSD::EndOfEvent(G4HCofThisEvent*)
{
}

//_____________________________________________________________________________
void MPPCSD::InitializeQESpline()
{
  // Use PDE data from KVC_OpticalProperties.hh instead of hardcoding a redundant copy here.
  auto graph = new TGraph(KVC_Optical::E_MPPC_PDE.size(), &KVC_Optical::E_MPPC_PDE[0], &KVC_Optical::R_MPPC_PDE[0]);
  m_qe_spline = new TSpline3("qe_spline", graph);
  m_range_min = m_qe_spline->GetXmin();
  m_range_max = m_qe_spline->GetXmax();

  delete graph;
}
