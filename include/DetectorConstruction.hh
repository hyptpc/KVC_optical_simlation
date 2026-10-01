// -*- C++ -*-

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <map>
#include <vector>

#include <G4String.hh>
#include <G4Types.hh>
#include <G4VUserDetectorConstruction.hh>

class DetectorMessenger; // Placeholder for a future messenger (not implemented yet)
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

//_____________________________________________________________________________
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

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
  std::vector<G4VPhysicalVolume*> m_gap_pvs;  // Air gap slabs (only if air_layer_thickness > 0)
  G4bool                          m_check_overlaps;

private:
  G4VPhysicalVolume* Construct() override;
  void ConstructElements();
  void ConstructMaterials();
  void ConstructKVC();
  void AddOpticalProperties();
  void AddSurfaceProperties();
  void DumpMaterialProperties(G4Material* mat); // For debugging (enabled with DEBUG)

  void CheckOverlaps(G4bool is_enabled) { m_check_overlaps = is_enabled; }
};

#endif
