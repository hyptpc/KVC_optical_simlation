#include "MPPCHit.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

G4Allocator<MPPCHit> MPPCHitAllocator;

MPPCHit::MPPCHit()
    : fPosition(G4ThreeVector()),
      fWorldPosition(G4ThreeVector()),
      fTime(0.),
      fEnergy(0.),
      fWaveLength(0.),
      fParticleID(0),
      fCopyNumber(0),
      fEventID(0),
      fDetectFlag(0)
{
}

MPPCHit::~MPPCHit() {}

MPPCHit::MPPCHit(const MPPCHit& right) : G4VHit() {
    fPosition = right.fPosition;
    fWorldPosition = right.fWorldPosition;
    fTime = right.fTime;
    fEnergy = right.fEnergy;
    fWaveLength = right.fWaveLength;
    fParticleID = right.fParticleID;
    fCopyNumber = right.fCopyNumber;
    fEventID = right.fEventID;
    fDetectFlag = right.fDetectFlag;
}

void MPPCHit::Print() const {
    G4cout << "MPPCHit: Position = " << fPosition
           << ", Time = " << G4BestUnit(fTime, "Time")
           << ", Energy = " << G4BestUnit(fEnergy, "Energy")
           << ", EventID = " << fEventID
           << G4endl;
}
