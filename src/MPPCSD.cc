// -*- C++ -*-

#include "MPPCSD.hh"

#include <G4EventManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4OpticalPhoton.hh>
#include <G4PhysicalConstants.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <Randomize.hh>

#include <TGraph.h>
#include <TSpline.h>

#include "ConfManager.hh"
#include "KVC_OpticalProperties.hh"

//_____________________________________________________________________________
MPPCSD::MPPCSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection(nullptr),
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
MPPCSD::~MPPCSD()
{
  delete m_qe_spline;
}

//_____________________________________________________________________________
void
MPPCSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<MPPCHit>(SensitiveDetectorName,
                                                     collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
MPPCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  // The step ends inside the MPPC: use the post-step point for the hit volume
  const auto post_step_point = aStep->GetPostStepPoint();
  const auto track = aStep->GetTrack();
  const auto definition = track->GetDefinition();
  const G4int particle_id = definition->GetPDGEncoding();
  if (definition != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

  const G4ThreeVector world_pos = post_step_point->GetPosition();
  const G4ThreeVector local_pos = post_step_point->GetTouchable()->GetHistory()
    ->GetTopTransform().TransformPoint(world_pos);
  const G4double hit_time = post_step_point->GetGlobalTime();
  const G4double energy = track->GetTotalEnergy();
  const G4double wave_length = (CLHEP::h_Planck * CLHEP::c_light / energy) / CLHEP::nm;
  const G4int copy_number = post_step_point->GetTouchableHandle()->GetCopyNumber();
  const G4int event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // -- kill track -----
  // Optical photons entering the MPPC are absorbed here regardless of the PDE result.
  track->SetTrackStatus(fStopAndKill);

  // -- PDE check -----
  G4double eval_energy = energy;
  if      (eval_energy < m_range_min) eval_energy = m_range_min;
  else if (eval_energy > m_range_max) eval_energy = m_range_max;

  G4double detection_prob = m_qe_spline->Eval(eval_energy) * m_qe_scale;
  if (detection_prob > 1.0) detection_prob = 1.0;

  const G4bool is_detected = (G4UniformRand() <= detection_prob);

  // -- record -----
  if (is_detected) {
    auto hit = new MPPCHit();
    hit->SetPosition(local_pos);
    hit->SetWorldPosition(world_pos);
    hit->SetEnergy(energy);
    hit->SetWaveLength(wave_length);
    hit->SetTime(hit_time);
    hit->SetParticleID(particle_id);
    hit->SetCopyNumber(copy_number);
    hit->SetEventID(event_id);
    hit->SetDetectFlag(1);

    m_hits_collection->insert(hit);
  }

  return true;
}

//_____________________________________________________________________________
void
MPPCSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
MPPCSD::InitializeQESpline()
{
  // Use the PDE data in KVC_OpticalProperties.hh
  auto graph = new TGraph(KVC_Optical::E_MPPC_PDE.size(),
                          &KVC_Optical::E_MPPC_PDE[0], &KVC_Optical::R_MPPC_PDE[0]);
  m_qe_spline = new TSpline3("qe_spline", graph);
  m_range_min = m_qe_spline->GetXmin();
  m_range_max = m_qe_spline->GetXmax();

  delete graph;
}
