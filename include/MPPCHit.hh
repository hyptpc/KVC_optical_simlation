// -*- C++ -*-

#ifndef MPPC_HIT_HH
#define MPPC_HIT_HH

#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

//_____________________________________________________________________________
class MPPCHit : public G4VHit
{
public:
  MPPCHit();
  ~MPPCHit() override;
  MPPCHit(const MPPCHit& right);

public:
  // Hit position (local coordinates of the MPPC)
  void          SetPosition(const G4ThreeVector& pos) { m_position = pos; }
  G4ThreeVector GetPosition() const { return m_position; }

  // Hit position (world coordinates)
  void          SetWorldPosition(const G4ThreeVector& pos) { m_world_position = pos; }
  G4ThreeVector GetWorldPosition() const { return m_world_position; }

  // Global time of the hit
  void     SetTime(G4double time) { m_time = time; }
  G4double GetTime() const { return m_time; }

  // Energy of the detected photon
  void     SetEnergy(G4double energy) { m_energy = energy; }
  G4double GetEnergy() const { return m_energy; }

  // Wavelength of the detected photon [nm]
  void     SetWaveLength(G4double wave_length) { m_wave_length = wave_length; }
  G4double GetWaveLength() const { return m_wave_length; }

  // PDG code of the particle
  void  SetParticleID(G4int pid) { m_particle_id = pid; }
  G4int GetParticleID() const { return m_particle_id; }

  // MPPC copy number (which MPPC was hit)
  void  SetCopyNumber(G4int copy_number) { m_copy_number = copy_number; }
  G4int GetCopyNumber() const { return m_copy_number; }

  // Event ID
  void  SetEventID(G4int event_id) { m_event_id = event_id; }
  G4int GetEventID() const { return m_event_id; }

  // Detection flag (1: detected)
  void  SetDetectFlag(G4int detect_flag) { m_detect_flag = detect_flag; }
  G4int GetDetectFlag() const { return m_detect_flag; }

  void Print() override;

private:
  G4ThreeVector m_position;       // Local position of the hit
  G4ThreeVector m_world_position; // World position of the hit
  G4double      m_time;           // Global time
  G4double      m_energy;         // Photon energy
  G4double      m_wave_length;    // Photon wavelength [nm]
  G4int         m_particle_id;    // PDG code
  G4int         m_copy_number;    // MPPC copy number
  G4int         m_event_id;       // Event ID
  G4int         m_detect_flag;    // Detection flag
};

#endif
