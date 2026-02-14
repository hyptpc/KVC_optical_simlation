#include "DetectorConstruction.hh"
#include "MPPCSD.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "KVC_OpticalProperties.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfManager.hh"

#define DEBUG 0

namespace
{
  auto& gConfMan = ConfManager::GetInstance();
}

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), m_check_overlaps(true),
    m_world_lv(nullptr), m_mother_lv(nullptr), m_blacksheet_lv(nullptr),
    m_mother_pv(nullptr), m_kvc_pv(nullptr), m_wrap_pv(nullptr)
{
}

//_____________________________________________________________________________
DetectorConstruction::~DetectorConstruction()
{
}

//_____________________________________________________________________________
G4VPhysicalVolume*
DetectorConstruction::Construct()
{
  using CLHEP::m;

  ConstructElements();
  ConstructMaterials();
  AddOpticalProperties();
  
  auto world_solid = new G4Box("WorldSolid", 1.*m/2, 1.*m/2, 1.*m/2);
  m_world_lv = new G4LogicalVolume(world_solid, m_material_map["Air"],
                                   "World");
  m_world_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto world_pv = new G4PVPlacement(nullptr, G4ThreeVector(), m_world_lv,
                                    "World", nullptr, false, 0, m_check_overlaps);

  ConstructKVC();
  AddSurfaceProperties();
  
  return world_pv;
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructElements()
{
  using CLHEP::g;
  using CLHEP::mole;
  /* G4Element(name, symbol, Z, A) */
  G4String name, symbol;
  G4double Z, A;
  name = "Hydrogen";
  m_element_map[name] = new G4Element(name, symbol="H",  Z=1.,
                                      A=1.00794 *g/mole);
  name = "Carbon";
  m_element_map[name] = new G4Element(name, symbol="C",  Z=6.,
                                      A=12.011 *g/mole);
  name = "Nitrogen";
  m_element_map[name] = new G4Element(name, symbol="N",  Z=7.,
                                      A=14.00674 *g/mole);
  name = "Oxygen";
  m_element_map[name] = new G4Element(name, symbol="O",  Z=8.,
                                      A=15.9994 *g/mole);
  name = "Sodium";
  m_element_map[name] = new G4Element(name, symbol="Na", Z=11.,
                                      A=22.989768 *g/mole);
  name = "Silicon";
  m_element_map[name] = new G4Element(name, symbol="Si", Z=14.,
                                      A=28.0855 *g/mole);
  name = "Phosphorus";
  m_element_map[name] = new G4Element(name, symbol="P", Z=15.,
                                      A=30.973762 *g/mole);
  name = "Sulfur";
  m_element_map[name] = new G4Element(name, symbol="S", Z=16.,
                                      A=32.066 *g/mole);
  name = "Chlorine";
  m_element_map[name] = new G4Element(name, symbol="Cl", Z=17.,
                                      A=35.453 *g/mole);
  name = "Argon";
  m_element_map[name] = new G4Element(name, symbol="Ar", Z=18.,
                                      A=39.948 *g/mole);
  name = "Potassium";
  m_element_map[name] = new G4Element(name, symbol="K", Z=19.,
                                      A=39.093 *g/mole);
  name = "Fluorine";
  m_element_map[name] = new G4Element(name, symbol="F", Z=9.,
                                      A=18.998 *g/mole);

  name = "Titanium";
  m_element_map[name] = new G4Element(name, symbol="Ti", Z=22.,
                                      A=47.867 * g/mole);

}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructMaterials()
{
  using CLHEP::g;
  using CLHEP::mg;
  using CLHEP::cm3;
  using CLHEP::mole;
  using CLHEP::STP_Temperature;

  /*
    G4Material(name, density, nelement, state, temperature, pressure);
    G4Material(name, z, a, density, state, temperature, pressure);
  */
  G4String name;
  G4double Z, A, density, massfraction;
  G4int natoms, nel, ncomponents;
  const G4double room_temp = STP_Temperature + 20.*CLHEP::kelvin;

  // Vacuum
  name = "Vacuum";
  m_material_map[name] =
    new G4Material(name, density=CLHEP::universe_mean_density, nel=2);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"], 0.7);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 0.3);
  
  // Air
  name = "Air";
  m_material_map[name] = new G4Material(name, density=1.2929e-03*g/cm3,
                                        nel=3, kStateGas, room_temp);
  G4double fracN  = 75.47;
  G4double fracO  = 23.20;
  G4double fracAr =  1.28;
  G4double denominator = fracN + fracO + fracAr;
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction=fracN/denominator);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction=fracO/denominator);
  m_material_map[name]->AddElement(m_element_map["Argon"],
                                   massfraction=fracAr/denominator);
  // Water
  name = "Water";
  m_material_map[name] = new G4Material(name, density=1.*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 1);

  // G10 epoxy glass
  name = "G10";
  m_material_map[name] = new G4Material(name, density=1.700*g/cm3,
                                        ncomponents=4);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"] , natoms=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"] , natoms=3);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"] , natoms=3);

  // Kapton
  name = "Kapton";
  m_material_map[name] = new G4Material(name, density=1.42*g/cm3,
                                        ncomponents=4);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],
                                   massfraction=0.0273);
  m_material_map[name]->AddElement(m_element_map["Carbon"],
                                   massfraction=0.7213);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction=0.0765);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction=0.1749);

  // Scintillator (Polystyene(C6H5CH=CH2))
  name = "Scintillator";
  m_material_map[name] = new G4Material(name, density=1.032*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=8);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=8);

  // EJ-232 (Plastic Scintillator, Polyvinyltoluene)
  name = "EJ232";
  m_material_map[name] = new G4Material(name, density=1.023*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=9);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=12);

  // CH2 Polyethelene
  name = "CH2";
  m_material_map[name] = new G4Material(name, density=0.95*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=2);

  // Aerogel
  name = "Aerogel";
  m_material_map[name] = new G4Material(name, density=0.2 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);
  
  // Quartz for KVC (SiO2, crystalline)
  name = "QuartzKVC";
  m_material_map[name] = new G4Material(name, density=2.64 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);

  // Acrylic for WC
  name = "Acrylic";
  m_material_map[name] = new G4Material(name, density=1.18 *g/cm3, nel=3);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=5);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],  natoms=8);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);

  // blacksheet
  name = "Blacksheet";
  m_material_map[name] = new G4Material(name, density=0.95 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],  natoms=2);

  // MPPC
  name = "MPPC";
  m_material_map[name] = new G4Material(name, 2.2 * g/cm3, 2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], 1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 2);

  // Epoxi
  name = "Epoxi";
  m_material_map[name] = new G4Material(name, density=1.1 *g/cm3, nel=4);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=21);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],  natoms=25);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=5);
  m_material_map[name]->AddElement(m_element_map["Chlorine"],  natoms=1);

  // Teflon
  name = "Teflon";
  m_material_map[name] = new G4Material(name, density=2.2 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=2);
  m_material_map[name]->AddElement(m_element_map["Fluorine"],  natoms=4);

  // Mylar
  name = "Mylar";
  m_material_map[name] = new G4Material(name, density=1.39*g/cm3, ncomponents=3);
  m_material_map[name]->AddElement(m_element_map["Carbon"], 10);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], 8);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 4);

  // EJ-510 White Reflective Paint (approximate)
  name = "EJ510";
  // Density: Representative value for white epoxy paint (literature/experience)
  m_material_map[name] = new G4Material(name, density = 1.6 * g/cm3, nel = 4);
  // Chemical composition (Epoxy + White pigment approximation)
  // Strict accuracy is not required (optics is dominated by surface properties)
  m_material_map[name]->AddElement(m_element_map["Carbon"],   natoms = 15);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms = 18);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],   natoms = 4);
  m_material_map[name]->AddElement(m_element_map["Titanium"], natoms = 1); // Representative of TiO2 pigment
}

//_____________________________________________________________________________
void
DetectorConstruction::AddOpticalProperties()
{
  using CLHEP::eV;
  using CLHEP::m;
  using CLHEP::mm;
  using CLHEP::cm;
  
  // +-----------------+
  // | Quartz Property |
  // +-----------------+
  auto quartz_prop = new G4MaterialPropertiesTable();
  quartz_prop->AddProperty("RINDEX", KVC_Optical::E_Quartz_RINDEX, KVC_Optical::R_Quartz_RINDEX);
  quartz_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Quartz_ABS, KVC_Optical::R_Quartz_ABS);
  m_material_map["QuartzKVC"]->SetMaterialPropertiesTable(quartz_prop);
  
  // +--------------+
  // | Air Property |
  // +--------------+
  auto air_prop = new G4MaterialPropertiesTable();
  air_prop->AddProperty("RINDEX", KVC_Optical::E_Air, KVC_Optical::R_Air_RINDEX);
  m_material_map["Air"]->SetMaterialPropertiesTable(air_prop);
  
  // +----------------------+
  // | Black sheet Property |
  // +----------------------+
  auto blacksheet_prop = new G4MaterialPropertiesTable();
  blacksheet_prop->AddProperty("RINDEX", KVC_Optical::E_Blacksheet, KVC_Optical::R_Blacksheet_RINDEX);
  blacksheet_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Blacksheet, KVC_Optical::R_Blacksheet_ABS);
  m_material_map["Blacksheet"]->SetMaterialPropertiesTable(blacksheet_prop);

  // +-----------------+
  // | Teflon Property |
  // +-----------------+
  auto teflon_prop = new G4MaterialPropertiesTable();
  teflon_prop->AddProperty("RINDEX", KVC_Optical::E_Teflon, KVC_Optical::R_Teflon_RINDEX);
  teflon_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Teflon, KVC_Optical::R_Teflon_ABS);
  m_material_map["Teflon"]->SetMaterialPropertiesTable(teflon_prop);

  // +----------------+
  // | Mylar Property |
  // +----------------+
  // Mylar surface is defined as dielectric_metal, so light does not penetrate. 
  // RINDEX and ABSLENGTH are defined here for potential future model updates.
  auto mylar_prop = new G4MaterialPropertiesTable();
  mylar_prop->AddProperty("RINDEX", KVC_Optical::E_Mylar, KVC_Optical::R_Mylar_RINDEX);
  mylar_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Mylar, KVC_Optical::R_Mylar_ABS);
  m_material_map["Mylar"]->SetMaterialPropertiesTable(mylar_prop);

  // +-----------------+
  // | EJ-510 Property |
  // +-----------------+
  // EJ-510 Property: Using estimated values for reflectivity grid.
  auto ej510_prop = new G4MaterialPropertiesTable();
  ej510_prop->AddProperty("RINDEX", KVC_Optical::E_EJ510_Bulk, KVC_Optical::R_EJ510_RINDEX);
  ej510_prop->AddProperty("ABSLENGTH", KVC_Optical::E_EJ510_Bulk, KVC_Optical::R_EJ510_ABS);
  m_material_map["EJ510"]->SetMaterialPropertiesTable(ej510_prop);
  
  // +---------------+
  // | MPPC Property |
  // +---------------+
  auto mppc_prop = new G4MaterialPropertiesTable();
  mppc_prop->AddProperty("RINDEX", KVC_Optical::E_MPPC, KVC_Optical::R_MPPC_RINDEX);
  // mppc_prop->AddProperty("ABSLENGTH", KVC_Optical::E_MPPC, KVC_Optical::R_MPPC_ABS);
  m_material_map["MPPC"]->SetMaterialPropertiesTable(mppc_prop);
  
  // +-------------------------------+
  // | MPPC surface (Epoxi) Property |
  // +-------------------------------+
  auto epoxi_prop = new G4MaterialPropertiesTable();
  epoxi_prop->AddProperty("RINDEX", KVC_Optical::E_Epoxi, KVC_Optical::R_Epoxi_RINDEX);
  epoxi_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Epoxi, KVC_Optical::R_Epoxi_ABS);
  m_material_map["Epoxi"]->SetMaterialPropertiesTable(epoxi_prop); 
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructKVC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::eV;

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();

  // Parameters from ConfManager 
  G4double quartz_thickness    = gConfMan.GetDouble("quartz_thickness") * mm;
  G4double air_layer_thickness = gConfMan.GetDouble("air_layer_thickness") * mm;
  G4double wrapper_thickness   = gConfMan.GetDouble("wrapper_thickness") * mm;
  G4int do_segmentize          = gConfMan.GetInt("do_segmentize");
  G4int wrap_type              = gConfMan.GetInt("wrap_type");

  G4ThreeVector kvc_size = (do_segmentize == 1)
    ? G4ThreeVector(26.0 * mm, 120.0 * mm, quartz_thickness)
    : G4ThreeVector(104.0 * mm, 120.0 * mm, quartz_thickness);

  G4ThreeVector origin_pos(0.0*mm, 0.0*mm, 0.0*mm);

  // Mother Volume (Air)
  auto mother_solid = new G4Box("KvcMotherSolid", 
                                kvc_size.x()/2.0 + 50.0*mm,
                                kvc_size.y()/2.0 + 50.0*mm,
                                kvc_size.z()/2.0 + 50.0*mm); 
  m_mother_lv = new G4LogicalVolume(mother_solid, m_material_map["Air"], "KvcMotherLV");
  m_mother_pv = new G4PVPlacement(nullptr, origin_pos, m_mother_lv,
                                  "KvcMotherPV", m_world_lv, false, 0, m_check_overlaps);
  m_mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Radiator
  auto kvc_solid = new G4Box("KvcSolid", 
			     kvc_size.x()/2.0,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0);
  auto kvc_lv = new G4LogicalVolume(kvc_solid, m_material_map["QuartzKVC"], "KvcLV");  
  m_kvc_pv = new G4PVPlacement(nullptr, origin_pos, kvc_lv, "KvcPV",
                               m_mother_lv, false, 0, m_check_overlaps);
  kvc_lv->SetVisAttributes(G4Colour::Yellow());

  // Wrapper
  G4Material* wrap_material = nullptr;
  if      (wrap_type == 0) wrap_material = m_material_map["Teflon"];
  else if (wrap_type == 1) wrap_material = m_material_map["Mylar"];
  else if (wrap_type == 2) wrap_material = m_material_map["EJ510"];
  else {
    G4Exception("DetectorConstruction::ConstructKVC", "InvalidWrapType", FatalException, "wrap_type must be 0,1,2");
  }

  auto wrap_solid_full = new G4Box("WrapSolidFull",
			     kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness);
  auto wrap_solid_cut  = new G4Box("WrapSolidCut",
			     kvc_size.x()/2.0 + air_layer_thickness,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0 + air_layer_thickness);
  
  G4SubtractionSolid* wrap_solid = new G4SubtractionSolid("WrapSolid", wrap_solid_full, wrap_solid_cut, nullptr, origin_pos);
  auto wrap_lv = new G4LogicalVolume(wrap_solid, wrap_material, "WrapLV");
  m_wrap_pv = new G4PVPlacement(nullptr, origin_pos, wrap_lv, "WrapPV", m_mother_lv, false, 0, m_check_overlaps); 
  wrap_lv->SetVisAttributes(G4Colour::White());

  // MPPC 
  G4ThreeVector mppc_size(6.0*mm, 6.0*mm, 1.0*mm);
  auto mppc_solid = new G4Box("MppcSolid", mppc_size.x()/2.0, mppc_size.y()/2.0, mppc_size.z()/2.0);
  auto mppc_lv = new G4LogicalVolume(mppc_solid, m_material_map["Epoxi"], "MppcLV");
  
  auto rot = new G4RotationMatrix;
  rot->rotateX(90.0*deg);
  G4int n_mppc = (do_segmentize == 1) ? 4 : 16;
  G4double offset = 0.0 * mm;

  if (6.0*mm < quartz_thickness && quartz_thickness < 12.0*mm) {
    for(G4int i=0; i<n_mppc; ++i){
      G4ThreeVector pos_up(-(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset, 0.0*mm);
      G4ThreeVector pos_low(-(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset, 0.0*mm);
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_up,  mppc_lv, "MppcPV", m_mother_lv, false, i,          m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_low, mppc_lv, "MppcPV", m_mother_lv, false, i+n_mppc,   m_check_overlaps));
    }
  } else if (12.0*mm <= quartz_thickness) {
    for(G4int i=0; i<n_mppc; ++i){
      G4double z_offset = quartz_thickness/6.0 + 1.0*mm;
      G4ThreeVector pos_up1( -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset,  z_offset);
      G4ThreeVector pos_up2( -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset, -z_offset);
      G4ThreeVector pos_low1(-(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset,  z_offset);
      G4ThreeVector pos_low2(-(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i), -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset, -z_offset);
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_up1,  mppc_lv, "MppcPV", m_mother_lv, false, i,          m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_up2,  mppc_lv, "MppcPV", m_mother_lv, false, i+n_mppc,   m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_low1, mppc_lv, "MppcPV", m_mother_lv, false, i+2*n_mppc, m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(rot, pos_low2, mppc_lv, "MppcPV", m_mother_lv, false, i+3*n_mppc, m_check_overlaps));
    }
  } else {
    G4Exception("DetectorConstruction::ConstructKVC", "InvalidQuartzThickness", FatalException, "Quartz thickness too small.");
  }
  mppc_lv->SetVisAttributes(G4Colour::Blue());
  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  mppc_lv->SetSensitiveDetector(mppcSD);

  // Blacksheet
  auto blacksheet_solid_full = new G4Box("BlacksheetSolidFull",
                                         kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness + 4.0*mm,
                                         kvc_size.y()/2.0 + 5.0*mm,
                                         kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness + 4.0*mm);
  auto blacksheet_solid_cut  = new G4Box("BlacksheetSolidCut",
                                         kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness + 1.0*mm,
                                         kvc_size.y()/2.0 + 2.0*mm,
                                         kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness + 1.0*mm);
  auto blacksheet_solid = new G4SubtractionSolid("BlacksheetSolid", blacksheet_solid_full, blacksheet_solid_cut, nullptr, origin_pos);
  m_blacksheet_lv = new G4LogicalVolume(blacksheet_solid, m_material_map["Blacksheet"], "BlacksheetLV");
  new G4PVPlacement(nullptr, origin_pos, m_blacksheet_lv, "BlacksheetPV", m_mother_lv, false, 0, m_check_overlaps);
  m_blacksheet_lv->SetVisAttributes(G4Colour::Black());
}

//_____________________________________________________________________________
void
DetectorConstruction::AddSurfaceProperties()
{
  using CLHEP::mm;
  using CLHEP::eV;

  G4int wrap_type              = gConfMan.GetInt("wrap_type");
  G4double air_layer_thickness = gConfMan.GetDouble("air_layer_thickness") * mm;
  G4int quartz_finish          = gConfMan.GetInt("quartz_finish");  // 0:polished, 1:ground
  G4double sigma_alpha         = 0.0;
  if (gConfMan.Check("Quartz_A_Alpha") && quartz_finish == 0) {
      sigma_alpha = gConfMan.GetDouble("Quartz_A_Alpha");
  } else if (gConfMan.Check("Quartz_B_Alpha") && quartz_finish == 1) {
      sigma_alpha = gConfMan.GetDouble("Quartz_B_Alpha");
  } else if (gConfMan.Check("sigma_alpha")) {
      sigma_alpha = gConfMan.GetDouble("sigma_alpha");
  }

  // Quartz Surface 
  auto surface_quartz = new G4OpticalSurface("surface_quartz");
  surface_quartz->SetModel(unified);
  surface_quartz->SetType(dielectric_dielectric);
  if(quartz_finish == 1){
    surface_quartz->SetFinish(ground);
  } else {
    surface_quartz->SetFinish(polished);
  }
  surface_quartz->SetSigmaAlpha(sigma_alpha);

  auto quartz_prop = new G4MaterialPropertiesTable();
  std::vector<G4double> e_surface = KVC_Optical::E_Unified_Surface;
  quartz_prop->AddConstProperty("SPECULARLOBECONSTANT",  gConfMan.GetDouble("quartz_specularLobe"), true);
  quartz_prop->AddConstProperty("SPECULARSPIKECONSTANT", gConfMan.GetDouble("quartz_specularSpike"), true);
  quartz_prop->AddConstProperty("BACKSCATTERCONSTANT",  gConfMan.GetDouble("quartz_backScatter"), true);

  G4double q_boundary_r = gConfMan.GetDouble("quartz_boundary_reflectivity");
  if (q_boundary_r >= 0.0) {
      quartz_prop->AddProperty("REFLECTIVITY", e_surface, std::vector<G4double>{q_boundary_r, q_boundary_r});
  }

  surface_quartz->SetMaterialPropertiesTable(quartz_prop);

  // Wrapper Surface
  G4OpticalSurface* wrap_surface = nullptr;

  // Wrappers (Teflon, Mylar, EJ-510)
  // Shared logic for wrapper configuration to avoid code duplication and confusion
  // Specific material properties are selected based on wrap_type

  G4OpticalSurface* surface_wrapper = new G4OpticalSurface("surface_wrapper");
  surface_wrapper->SetModel(unified);

  auto wrapper_prop = new G4MaterialPropertiesTable();

  // Common parameters that apply to all wrappers (or at least the ones using unified model)
  // These keys are "wrapper_..." to indicate they are tunable parameters for whatever wrapper is selected.
  // Note: Mylar might ignore some of these due to being dielectric_metal / polished.
  
  if (wrap_type == 0) { // Teflon
      surface_wrapper->SetType(dielectric_dielectric);
      surface_wrapper->SetFinish(groundfrontpainted);
      surface_wrapper->SetSigmaAlpha(gConfMan.GetDouble("teflon_sigma_alpha")); 

      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_PTFE_Thin);
      wrapper_prop->AddConstProperty("SPECULARLOBECONSTANT",  gConfMan.GetDouble("teflon_specularLobe"), true);
      wrapper_prop->AddConstProperty("SPECULARSPIKECONSTANT", gConfMan.GetDouble("teflon_specularSpike"), true);
      wrapper_prop->AddConstProperty("BACKSCATTERCONSTANT",   gConfMan.GetDouble("teflon_backScatter"), true);

  } else if (wrap_type == 1) { // Mylar
      surface_wrapper->SetType(dielectric_metal);
      surface_wrapper->SetFinish(polished);
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_AlMylar);

  } else if (wrap_type == 2) { // EJ-510
      surface_wrapper->SetType(dielectric_dielectric);
      surface_wrapper->SetFinish(groundfrontpainted);
      surface_wrapper->SetSigmaAlpha(gConfMan.GetDouble("ej510_sigma_alpha")); 

      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_EJ510);
      wrapper_prop->AddConstProperty("SPECULARLOBECONSTANT",  gConfMan.GetDouble("ej510_specularLobe"), true);
      wrapper_prop->AddConstProperty("SPECULARSPIKECONSTANT", gConfMan.GetDouble("ej510_specularSpike"), true);
      wrapper_prop->AddConstProperty("BACKSCATTERCONSTANT",   gConfMan.GetDouble("ej510_backScatter"), true);

  } else {
       G4Exception("DetectorConstruction::AddSurfaceProperties", "InvalidWrapType", FatalException, "wrap_type must be 0,1,2");
  }

  surface_wrapper->SetMaterialPropertiesTable(wrapper_prop);

  // Assign the created wrapper surface
  wrap_surface = surface_wrapper;

  // Border Surfaces
  if (m_kvc_pv && m_mother_pv && m_wrap_pv) {
    if (air_layer_thickness > 0.0) {
      new G4LogicalBorderSurface("QuartzToAir", m_kvc_pv,    m_mother_pv, surface_quartz);
      new G4LogicalBorderSurface("AirToQuartz", m_mother_pv, m_kvc_pv,    surface_quartz);
      new G4LogicalBorderSurface("AirToWrap",   m_mother_pv, m_wrap_pv,   wrap_surface);
      new G4LogicalBorderSurface("WrapToAir",   m_wrap_pv,   m_mother_pv, wrap_surface);
    } else {
      new G4LogicalBorderSurface("QuartzToWrap", m_kvc_pv,  m_wrap_pv, wrap_surface);
      new G4LogicalBorderSurface("WrapToQuartz", m_wrap_pv, m_kvc_pv,  wrap_surface);
    }
  }

  // BlackSheet Surface
  auto surface_bs = new G4OpticalSurface("surface_bs", unified, ground, dielectric_metal);
  auto bs_prop = new G4MaterialPropertiesTable();
  bs_prop->AddProperty("REFLECTIVITY", KVC_Optical::E_Blacksheet, KVC_Optical::R_Blacksheet_REFLECTIVITY);
  surface_bs->SetMaterialPropertiesTable(bs_prop);
  if (m_blacksheet_lv) new G4LogicalSkinSurface("BlackSheetSurface", m_blacksheet_lv, surface_bs);

  // MPPC Surface (Added for Method A)
  // Retrieve MppcLV from Store since it's not a member
  auto mppc_lv = G4LogicalVolumeStore::GetInstance()->GetVolume("MppcLV", false);
  if (mppc_lv) {
      auto surface_mppc = new G4OpticalSurface("surface_mppc");
      surface_mppc->SetType(dielectric_dielectric);
      surface_mppc->SetFinish(polished);
      surface_mppc->SetModel(unified); // or glisur

      auto mppc_surf_prop = new G4MaterialPropertiesTable();
#ifdef USE_SURFACE_PDE
      // Method A: Add EFFICIENCY to the surface
      // Use the PDE data from KVC_OpticalProperties (assuming it's roughly the same as manufacturer PDE)
      // Note: If Manufacturer PDE includes surface reflection, and we simulate surface reflection here,
      // strictly speaking we should strip reflection from PDE?
      // Usually manufacturer PDE (e.g. 40%) is "Probability of detection per incident photon".
      // Geant4 EFFICIENCY is "Probability of detection IF absorbed at boundary".
      // If we set type=dielectric_dielectric, Fresnel reflection happens.
      // So some photons reflect. Absorbed ones get checked against EFFICIENCY.
      // So EFFICIENCY should be ~ PDE / (1 - R).
      // However, for simplicity and standard usage, often raw PDE is used or R is set to 0.
      // Here we use the raw PDE data as EFFICIENCY and let Fresnel handle reflection.
      // This might slightly underestimate if PDE was "per incident" (which includes reflection loss).
      // But user report says "Method A is recommended".
      
      
      mppc_surf_prop->AddProperty("EFFICIENCY", KVC_Optical::E_MPPC_PDE, KVC_Optical::R_MPPC_PDE);
#endif
      surface_mppc->SetMaterialPropertiesTable(mppc_surf_prop);
      new G4LogicalSkinSurface("MppcSurface", mppc_lv, surface_mppc);

      // Direct Optical Coupling (Border Surface between Quartz and MPPC)
      if (m_kvc_pv) {
          for (size_t i = 0; i < m_mppc_pvs.size(); ++i) {
              new G4LogicalBorderSurface("QuartzToMppc", m_kvc_pv, m_mppc_pvs[i], surface_mppc);
              new G4LogicalBorderSurface("MppcToQuartz", m_mppc_pvs[i], m_kvc_pv, surface_mppc);
          }
      }
  }
}

//_____________________________________________________________________________
void DetectorConstruction::DumpMaterialProperties(G4Material* mat)
{
#if DEBUG
  G4cout << "=== Material: " << mat->GetName() << " ===" << G4endl;

  auto matPropTable = mat->GetMaterialPropertiesTable();
  if (!matPropTable) {
    G4cout << "No material properties table found." << G4endl;
    return;
  }

  std::vector<G4String> propertyNames = {"RINDEX", "ABSLENGTH", "REFLECTIVITY"};

  for (const auto& prop : propertyNames) {
    if (matPropTable->ConstPropertyExists(prop)) {
      G4cout << prop << ": " << matPropTable->GetConstProperty(prop) << G4endl;
    }
  }

  for (const auto& prop : propertyNames) {
    if (matPropTable->GetProperty(prop)) {
      G4MaterialPropertyVector* mpv = matPropTable->GetProperty(prop);
      G4cout << prop << ":" << G4endl;
      for (size_t i = 0; i < mpv->GetVectorLength(); i++) {
        G4cout << "Energy: " << mpv->Energy(i) / eV << " eV, "
               << " Value: " << (*mpv)[i] << G4endl;
      }
    }
  }
#endif
}
