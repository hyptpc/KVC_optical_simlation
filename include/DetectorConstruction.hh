#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4Element.hh"
#include "G4Material.hh"

#include <vector>
#include "G4VPhysicalVolume.hh"

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

private:
  std::map<G4String, G4Element*>  m_element_map;
  std::map<G4String, G4Material*> m_material_map;
  G4LogicalVolume*                m_world_lv;
  G4LogicalVolume*                m_mother_lv;
  G4LogicalVolume*                m_blacksheet_lv;
  G4VPhysicalVolume*              m_mother_pv;
  G4VPhysicalVolume*              m_kvc_pv;
  G4VPhysicalVolume*              m_wrap_pv;
  std::vector<G4VPhysicalVolume*> m_mppc_pvs;
  G4bool                          m_check_overlaps;

private:
  virtual G4VPhysicalVolume* Construct();
  void ConstructElements();
  void ConstructMaterials();
  void ConstructKVC();
  void AddOpticalProperties();
  void AddSurfaceProperties();
  void DumpMaterialProperties(G4Material* mat);

  void CheckOverlaps(G4bool flag) { m_check_overlaps = flag; }
};

#endif /*OpNoviceDetectorConstruction_h*/
