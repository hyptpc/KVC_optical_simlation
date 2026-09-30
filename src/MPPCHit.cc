// -*- C++ -*-

#include "MPPCHit.hh"

#include <G4UnitsTable.hh>
#include <G4ios.hh>

//_____________________________________________________________________________
MPPCHit::MPPCHit()
  : G4VHit(),
    m_position(),
    m_world_position(),
    m_time(0.),
    m_energy(0.),
    m_wave_length(0.),
    m_particle_id(0),
    m_copy_number(0),
    m_event_id(0),
    m_detect_flag(0)
{
}

//_____________________________________________________________________________
MPPCHit::~MPPCHit()
{
}

//_____________________________________________________________________________
MPPCHit::MPPCHit(const MPPCHit& right)
  : G4VHit(),
    m_position(right.m_position),
    m_world_position(right.m_world_position),
    m_time(right.m_time),
    m_energy(right.m_energy),
    m_wave_length(right.m_wave_length),
    m_particle_id(right.m_particle_id),
    m_copy_number(right.m_copy_number),
    m_event_id(right.m_event_id),
    m_detect_flag(right.m_detect_flag)
{
}

//_____________________________________________________________________________
void
MPPCHit::Print()
{
  G4cout << "MPPCHit: Position = " << m_position
         << ", Time = " << G4BestUnit(m_time, "Time")
         << ", Energy = " << G4BestUnit(m_energy, "Energy")
         << ", EventID = " << m_event_id
         << G4endl;
}
