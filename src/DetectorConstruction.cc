// -*- C++ -*-

#include "DetectorConstruction.hh"

#include <cmath>
#include <sstream>
#include <vector>

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Element.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4VisAttributes.hh>

#include "ConfManager.hh"
#include "KVC_OpticalProperties.hh"
#include "MPPCSD.hh"

#define DEBUG 0

namespace
{
  auto& gConfMan = ConfManager::GetInstance();

  //___________________________________________________________________________
  // Register the unified-model constants of an optical surface from the conf keys
  //   <material>_specularSpike, <material>_specularLobe,
  //   <material>_backScatter, <material>_diffuseLobe
  // - The three probabilities are registered as (constant) vector properties,
  //   which is how G4OpBoundaryProcess reads them.
  // - The diffuse lobe is the remainder 1 - (spike + lobe + backscatter);
  //   DIFFUSELOBECONSTANT is not a Geant4 property and is not read.
  //   If the four values do not sum to 1, a warning is issued.
  void
  RegisterUnifiedConstants(G4MaterialPropertiesTable* prop, const G4String& material)
  {
    const G4double spike   = gConfMan.GetDouble(material + "_specularSpike");
    const G4double lobe    = gConfMan.GetDouble(material + "_specularLobe");
    const G4double back    = gConfMan.GetDouble(material + "_backScatter");
    const G4double diffuse = gConfMan.GetDouble(material + "_diffuseLobe");

    const auto& energy = KVC_Optical::E_Unified_Surface;
    const auto constant = [&energy](G4double value) {
      return std::vector<G4double>(energy.size(), value);
    };
    prop->AddProperty("SPECULARSPIKECONSTANT", energy, constant(spike));
    prop->AddProperty("SPECULARLOBECONSTANT", energy, constant(lobe));
    prop->AddProperty("BACKSCATTERCONSTANT", energy, constant(back));
    prop->AddConstProperty("DIFFUSELOBECONSTANT", diffuse, true);

    constexpr G4double tolerance = 1.e-6;
    const G4double sum = spike + lobe + back + diffuse;
    if (std::abs(sum - 1.0) > tolerance) {
      std::ostringstream message;
      message << "Sum of the " << material << " unified-model constants"
              << " (specularSpike + specularLobe + backScatter + diffuseLobe) is "
              << sum << ", not 1." << G4endl
              << "Geant4 uses the diffuse lobe = 1 - (specularSpike + specularLobe + backScatter) = "
              << 1.0 - (spike + lobe + back) << " (" << material << "_diffuseLobe = "
              << diffuse << " is not used).";
      G4Exception("DetectorConstruction::AddSurfaceProperties", "UnifiedConstantsSum",
                  JustWarning, message);
    }
  }

  //___________________________________________________________________________
  // Set the unified-model parameters of an optical surface from the conf keys
  //   <material>_sigma_alpha and the unified-model constants (see above)
  // in the same format for every surface. Depending on the surface type and
  // finish, Geant4 ignores some of them (see README).
  void
  SetUnifiedParameters(G4OpticalSurface* surface, G4MaterialPropertiesTable* prop,
                       const G4String& material)
  {
    surface->SetSigmaAlpha(gConfMan.GetDouble(material + "_sigma_alpha"));
    RegisterUnifiedConstants(prop, material);
  }
}

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    m_element_map(),
    m_material_map(),
    m_world_lv(nullptr),
    m_mother_lv(nullptr),
    m_blacksheet_lv(nullptr),
    m_mother_pv(nullptr),
    m_kvc_pv(nullptr),
    m_wrap_pv(nullptr),
    m_mppc_pvs(),
    m_gap_pvs(),
    m_check_overlaps(true)
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
  G4double density, massfraction;
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

  std::vector<G4double> r_quartz_abs = KVC_Optical::R_Quartz_ABS;
  if (gConfMan.Check("quartz_abs_scale")) {
    const G4double abs_scale = gConfMan.GetDouble("quartz_abs_scale");
    for (auto& abs_length : r_quartz_abs) abs_length *= abs_scale;
  }
  quartz_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Quartz_ABS, r_quartz_abs);
  m_material_map["QuartzKVC"]->SetMaterialPropertiesTable(quartz_prop);

  // +--------------+
  // | Air Property |
  // +--------------+
  G4double air_rindex = 1.0;
  if (gConfMan.Check("air_rindex")) air_rindex = gConfMan.GetDouble("air_rindex");
  auto air_prop = new G4MaterialPropertiesTable();
  air_prop->AddProperty("RINDEX", KVC_Optical::E_Air, std::vector<G4double>{air_rindex, air_rindex});
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
  const G4int wrap_type = gConfMan.GetInt("wrap_type");
  G4double teflon_rindex = 1.35;
  if (gConfMan.Check("teflon_rindex")) teflon_rindex = gConfMan.GetDouble("teflon_rindex");

  auto teflon_prop = new G4MaterialPropertiesTable();
  teflon_prop->AddProperty("RINDEX", KVC_Optical::E_Teflon, std::vector<G4double>{teflon_rindex, teflon_rindex});

  if (wrap_type == 3) {
    // Transmissive Teflon: long absorption length
    const std::vector<G4double> abs_long = { 10.0*m, 10.0*m };
    teflon_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Teflon, abs_long);
  } else {
    // Standard Teflon: opaque / absorptive bulk
    teflon_prop->AddProperty("ABSLENGTH", KVC_Optical::E_Teflon, KVC_Optical::R_Teflon_ABS);
  }
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
  using CLHEP::deg;
  using CLHEP::mm;

  // Parameters from ConfManager
  const G4double quartz_thickness    = gConfMan.GetDouble("quartz_thickness") * mm;
  const G4double air_layer_thickness = gConfMan.GetDouble("air_layer_thickness") * mm;
  const G4double wrapper_thickness   = gConfMan.GetDouble("wrapper_thickness") * mm;
  const G4int    do_segmentize       = gConfMan.GetInt("do_segmentize");
  const G4int    wrap_type           = gConfMan.GetInt("wrap_type");

  const G4ThreeVector kvc_size = (do_segmentize == 1)
    ? G4ThreeVector(26.0 * mm, 120.0 * mm, quartz_thickness)
    : G4ThreeVector(104.0 * mm, 120.0 * mm, quartz_thickness);

  const G4ThreeVector origin_pos(0.0*mm, 0.0*mm, 0.0*mm);

  // Mother volume (air)
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
  else if (wrap_type == 3) wrap_material = m_material_map["Teflon"]; // Transmissive Teflon
  else {
    G4Exception("DetectorConstruction::ConstructKVC", "InvalidWrapType",
                FatalException, "wrap_type must be 0,1,2,3");
  }

  auto wrap_solid_full = new G4Box("WrapSolidFull",
                                   kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness,
                                   kvc_size.y()/2.0,
                                   kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness);
  auto wrap_solid_cut  = new G4Box("WrapSolidCut",
                                   kvc_size.x()/2.0 + air_layer_thickness,
                                   kvc_size.y()/2.0,
                                   kvc_size.z()/2.0 + air_layer_thickness);
  auto wrap_solid = new G4SubtractionSolid("WrapSolid", wrap_solid_full, wrap_solid_cut,
                                           nullptr, origin_pos);
  auto wrap_lv = new G4LogicalVolume(wrap_solid, wrap_material, "WrapLV");
  m_wrap_pv = new G4PVPlacement(nullptr, origin_pos, wrap_lv, "WrapPV",
                                m_mother_lv, false, 0, m_check_overlaps);
  wrap_lv->SetVisAttributes(G4Colour::White());

  // Air gap between the quartz and the wrapper, made of four slabs (x and z sides), so
  // that the upper / lower (y) end faces of the quartz touch the mother volume directly
  // and can have a different optical surface (see AddSurfaceProperties).
  if (air_layer_thickness > 0.0) {
    auto gap_x_solid = new G4Box("AirGapXSolid", air_layer_thickness/2.0, kvc_size.y()/2.0,
                                 kvc_size.z()/2.0 + air_layer_thickness);
    auto gap_z_solid = new G4Box("AirGapZSolid", kvc_size.x()/2.0, kvc_size.y()/2.0,
                                 air_layer_thickness/2.0);
    auto gap_x_lv = new G4LogicalVolume(gap_x_solid, m_material_map["Air"], "AirGapXLV");
    auto gap_z_lv = new G4LogicalVolume(gap_z_solid, m_material_map["Air"], "AirGapZLV");
    gap_x_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
    gap_z_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
    const G4double gap_x_pos = kvc_size.x()/2.0 + air_layer_thickness/2.0;
    const G4double gap_z_pos = kvc_size.z()/2.0 + air_layer_thickness/2.0;
    m_gap_pvs.push_back(new G4PVPlacement(nullptr, G4ThreeVector( gap_x_pos, 0., 0.), gap_x_lv, "AirGapPV", m_mother_lv, false, 0, m_check_overlaps));
    m_gap_pvs.push_back(new G4PVPlacement(nullptr, G4ThreeVector(-gap_x_pos, 0., 0.), gap_x_lv, "AirGapPV", m_mother_lv, false, 1, m_check_overlaps));
    m_gap_pvs.push_back(new G4PVPlacement(nullptr, G4ThreeVector(0., 0.,  gap_z_pos), gap_z_lv, "AirGapPV", m_mother_lv, false, 2, m_check_overlaps));
    m_gap_pvs.push_back(new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -gap_z_pos), gap_z_lv, "AirGapPV", m_mother_lv, false, 3, m_check_overlaps));
  }

  // MPPC
  const G4ThreeVector mppc_size(6.0*mm, 6.0*mm, 1.0*mm);
  auto mppc_solid = new G4Box("MppcSolid", mppc_size.x()/2.0, mppc_size.y()/2.0, mppc_size.z()/2.0);
  auto mppc_lv = new G4LogicalVolume(mppc_solid, m_material_map["Epoxi"], "MppcLV");

  auto mppc_rot = new G4RotationMatrix;
  mppc_rot->rotateX(90.0*deg);
  const G4int n_mppc = (do_segmentize == 1) ? 4 : 16; // MPPCs per row
  const G4double mppc_gap = 0.0 * mm;                 // Gap between the quartz and the MPPC
  const G4double mppc_pitch = mppc_size.x() + 0.5*mm;
  const G4double y_up  =  kvc_size.y()/2.0 + mppc_size.z()/2.0 + mppc_gap;
  const G4double y_low = -kvc_size.y()/2.0 - mppc_size.z()/2.0 - mppc_gap;

  if (6.0*mm < quartz_thickness && quartz_thickness < 12.0*mm) {
    // One row of MPPCs on each of the upper and lower faces
    for (G4int i = 0; i < n_mppc; ++i) {
      const G4double x = -mppc_pitch * ((n_mppc-1)/2.0 - i);
      const G4ThreeVector pos_up(x, y_up, 0.0*mm);
      const G4ThreeVector pos_low(x, y_low, 0.0*mm);
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_up,  mppc_lv, "MppcPV", m_mother_lv, false, i,          m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_low, mppc_lv, "MppcPV", m_mother_lv, false, i+n_mppc,   m_check_overlaps));
    }
  } else if (12.0*mm <= quartz_thickness) {
    // Two rows of MPPCs on each of the upper and lower faces
    const G4double z_offset = quartz_thickness/6.0 + 1.0*mm;
    for (G4int i = 0; i < n_mppc; ++i) {
      const G4double x = -mppc_pitch * ((n_mppc-1)/2.0 - i);
      const G4ThreeVector pos_up1( x, y_up,   z_offset);
      const G4ThreeVector pos_up2( x, y_up,  -z_offset);
      const G4ThreeVector pos_low1(x, y_low,  z_offset);
      const G4ThreeVector pos_low2(x, y_low, -z_offset);
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_up1,  mppc_lv, "MppcPV", m_mother_lv, false, i,          m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_up2,  mppc_lv, "MppcPV", m_mother_lv, false, i+n_mppc,   m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_low1, mppc_lv, "MppcPV", m_mother_lv, false, i+2*n_mppc, m_check_overlaps));
      m_mppc_pvs.push_back(new G4PVPlacement(mppc_rot, pos_low2, mppc_lv, "MppcPV", m_mother_lv, false, i+3*n_mppc, m_check_overlaps));
    }
  } else {
    G4Exception("DetectorConstruction::ConstructKVC", "InvalidQuartzThickness",
                FatalException, "Quartz thickness too small.");
  }
  mppc_lv->SetVisAttributes(G4Colour::Blue());
  auto mppc_sd = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppc_sd);
  mppc_lv->SetSensitiveDetector(mppc_sd);

  // Blacksheet
  auto blacksheet_solid_full = new G4Box("BlacksheetSolidFull",
                                         kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness + 4.0*mm,
                                         kvc_size.y()/2.0 + 5.0*mm,
                                         kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness + 4.0*mm);
  auto blacksheet_solid_cut  = new G4Box("BlacksheetSolidCut",
                                         kvc_size.x()/2.0 + air_layer_thickness + wrapper_thickness + 1.0*mm,
                                         kvc_size.y()/2.0 + 2.0*mm,
                                         kvc_size.z()/2.0 + air_layer_thickness + wrapper_thickness + 1.0*mm);
  auto blacksheet_solid = new G4SubtractionSolid("BlacksheetSolid", blacksheet_solid_full,
                                                 blacksheet_solid_cut, nullptr, origin_pos);
  m_blacksheet_lv = new G4LogicalVolume(blacksheet_solid, m_material_map["Blacksheet"], "BlacksheetLV");
  new G4PVPlacement(nullptr, origin_pos, m_blacksheet_lv, "BlacksheetPV",
                    m_mother_lv, false, 0, m_check_overlaps);
  m_blacksheet_lv->SetVisAttributes(G4Colour::Black());
}

//_____________________________________________________________________________
void
DetectorConstruction::AddSurfaceProperties()
{
  using CLHEP::mm;

  const G4int    wrap_type           = gConfMan.GetInt("wrap_type");
  const G4double air_layer_thickness = gConfMan.GetDouble("air_layer_thickness") * mm;
  const G4int    quartz_finish       = gConfMan.GetInt("quartz_finish"); // 0: polished, 1: ground
  G4double quartz_sigma_alpha = 0.0;
  if (gConfMan.Check("Quartz_A_Alpha") && quartz_finish == 0) {
    quartz_sigma_alpha = gConfMan.GetDouble("Quartz_A_Alpha");
  } else if (gConfMan.Check("Quartz_B_Alpha") && quartz_finish == 1) {
    quartz_sigma_alpha = gConfMan.GetDouble("Quartz_B_Alpha");
  } else if (gConfMan.Check("sigma_alpha")) {
    quartz_sigma_alpha = gConfMan.GetDouble("sigma_alpha");
  }

  // Quartz surface (used ONLY for the quartz-air gap boundaries;
  // the Quartz-MPPC boundary uses surface_mppc_refl below, i.e. polished / mirror-like)
  auto surface_quartz = new G4OpticalSurface("surface_quartz");
  surface_quartz->SetModel(unified);
  surface_quartz->SetType(dielectric_dielectric);
  if (quartz_finish == 1) {
    surface_quartz->SetFinish(ground);
  } else {
    surface_quartz->SetFinish(polished);
  }

  auto quartz_prop = new G4MaterialPropertiesTable();
  const std::vector<G4double> e_surface = KVC_Optical::E_Unified_Surface;
  // The quartz sigma_alpha is selected by quartz_finish (Quartz_A_Alpha / Quartz_B_Alpha)
  surface_quartz->SetSigmaAlpha(quartz_sigma_alpha);
  RegisterUnifiedConstants(quartz_prop, "quartz");

  const G4double quartz_reflectivity = gConfMan.GetDouble("quartz_boundary_reflectivity");
  if (quartz_reflectivity >= 0.0) {
    quartz_prop->AddProperty("REFLECTIVITY", e_surface,
                             std::vector<G4double>{quartz_reflectivity, quartz_reflectivity});
  }

  surface_quartz->SetMaterialPropertiesTable(quartz_prop);

  // Upper / lower (y) end faces of the quartz (MPPC side): always polished, also for the
  // frosted quartz B
  auto surface_quartz_end = new G4OpticalSurface("surface_quartz_end");
  surface_quartz_end->SetModel(unified);
  surface_quartz_end->SetType(dielectric_dielectric);
  surface_quartz_end->SetFinish(polished);
  auto quartz_end_prop = new G4MaterialPropertiesTable();
  if (quartz_reflectivity >= 0.0) {
    quartz_end_prop->AddProperty("REFLECTIVITY", e_surface,
                                 std::vector<G4double>{quartz_reflectivity, quartz_reflectivity});
  }
  surface_quartz_end->SetMaterialPropertiesTable(quartz_end_prop);

  // Wrapper surface (Teflon, Mylar, EJ-510), selected by wrap_type
  auto surface_wrapper = new G4OpticalSurface("surface_wrapper");
  surface_wrapper->SetModel(unified);
  auto wrapper_prop = new G4MaterialPropertiesTable();

  // Material whose conf keys (<material>_sigma_alpha, ...) are used for the wrapper surface
  const G4bool is_teflon = gConfMan.Check("is_teflon") && gConfMan.GetInt("is_teflon") == 1;
  const G4bool is_paint  = gConfMan.Check("is_paint")  && gConfMan.GetInt("is_paint")  == 1;

  if (wrap_type == 0) { // Teflon
    surface_wrapper->SetType(dielectric_dielectric);
    surface_wrapper->SetFinish(groundfrontpainted);
    SetUnifiedParameters(surface_wrapper, wrapper_prop, "teflon");

    std::vector<G4double> r_ptfe = KVC_Optical::R_PTFE_Thin;
    const G4double r_scale = gConfMan.GetDouble("teflon_reflectivity_scale");
    for (auto& r : r_ptfe) r *= r_scale;
    wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, r_ptfe);

  } else if (wrap_type == 1) { // Specular wrapper (Mylar, Teflon, or paint)
    surface_wrapper->SetType(dielectric_metal);

    if (is_teflon) {
      surface_wrapper->SetFinish(ground);
      SetUnifiedParameters(surface_wrapper, wrapper_prop, "teflon");

      const G4double r_scale = gConfMan.GetDouble("teflon_reflectivity_scale");
      std::vector<G4double> r_ptfe = KVC_Optical::R_PTFE_Thin;
      for (auto& r : r_ptfe) r *= r_scale;
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, r_ptfe);
    } else if (is_paint) {
      surface_wrapper->SetFinish(ground);
      SetUnifiedParameters(surface_wrapper, wrapper_prop, "ej510");
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_EJ510);
    } else {
      // Default: aluminized Mylar (polished mirror, no tunable surface parameters)
      surface_wrapper->SetFinish(polished);
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_AlMylar);
    }

  } else if (wrap_type == 2) { // EJ-510 style volume reflection (paint, or Teflon)
    surface_wrapper->SetType(dielectric_dielectric);
    surface_wrapper->SetFinish(groundfrontpainted);

    if (is_teflon) {
      SetUnifiedParameters(surface_wrapper, wrapper_prop, "teflon");

      const G4double r_scale = gConfMan.GetDouble("teflon_reflectivity_scale");
      std::vector<G4double> r_vec = KVC_Optical::R_EJ510; // Use the paint grid as the base
      for (auto& r : r_vec) r *= r_scale;
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, r_vec);
    } else {
      SetUnifiedParameters(surface_wrapper, wrapper_prop, "ej510");
      wrapper_prop->AddProperty("REFLECTIVITY", KVC_Optical::Energy, KVC_Optical::R_EJ510);
    }

  } else if (wrap_type == 3) { // Transmissive Teflon
    // Model: dielectric_dielectric + ground (rough interface).
    // Light can enter the Teflon volume based on Fresnel / micro-facets.
    // Requires the Teflon RINDEX and a long ABSLENGTH (set in AddOpticalProperties).
    // Do NOT set REFLECTIVITY here: Fresnel handles reflection vs transmission.
    surface_wrapper->SetType(dielectric_dielectric);
    surface_wrapper->SetFinish(ground);
    SetUnifiedParameters(surface_wrapper, wrapper_prop, "teflon");

  } else {
    G4Exception("DetectorConstruction::AddSurfaceProperties", "InvalidWrapType",
                FatalException, "wrap_type must be 0,1,2,3");
  }

  surface_wrapper->SetMaterialPropertiesTable(wrapper_prop);

  // Border surfaces
  if (m_kvc_pv && m_mother_pv && m_wrap_pv) {
    // End faces of the quartz (outside the MPPCs) <-> air: polished
    new G4LogicalBorderSurface("QuartzToAir", m_kvc_pv,    m_mother_pv, surface_quartz_end);
    new G4LogicalBorderSurface("AirToQuartz", m_mother_pv, m_kvc_pv,    surface_quartz_end);
    if (air_layer_thickness > 0.0) {
      for (const auto gap_pv : m_gap_pvs) {
        new G4LogicalBorderSurface("QuartzToGap", m_kvc_pv, gap_pv,    surface_quartz);
        new G4LogicalBorderSurface("GapToQuartz", gap_pv,   m_kvc_pv,  surface_quartz);
        new G4LogicalBorderSurface("GapToWrap",   gap_pv,   m_wrap_pv, surface_wrapper);
        new G4LogicalBorderSurface("WrapToGap",   m_wrap_pv, gap_pv,   surface_wrapper);
      }
      new G4LogicalBorderSurface("AirToWrap",   m_mother_pv, m_wrap_pv,   surface_wrapper);
      new G4LogicalBorderSurface("WrapToAir",   m_wrap_pv,   m_mother_pv, surface_wrapper);
    } else {
      new G4LogicalBorderSurface("QuartzToWrap", m_kvc_pv,  m_wrap_pv, surface_wrapper);
      new G4LogicalBorderSurface("WrapToQuartz", m_wrap_pv, m_kvc_pv,  surface_wrapper);
    }
  }

  // Blacksheet surface
  auto surface_blacksheet = new G4OpticalSurface("surface_bs", unified, ground, dielectric_metal);
  auto blacksheet_prop = new G4MaterialPropertiesTable();
  blacksheet_prop->AddProperty("REFLECTIVITY", KVC_Optical::E_Blacksheet,
                               KVC_Optical::R_Blacksheet_REFLECTIVITY);
  surface_blacksheet->SetMaterialPropertiesTable(blacksheet_prop);
  if (m_blacksheet_lv) {
    new G4LogicalSkinSurface("BlackSheetSurface", m_blacksheet_lv, surface_blacksheet);
  }

  // Optical boundary between the quartz and the MPPCs for Fresnel reflection.
  // Photon detection itself is handled in MPPCSD.
  auto mppc_lv = G4LogicalVolumeStore::GetInstance()->GetVolume("MppcLV", false);
  if (mppc_lv) {
    auto surface_mppc_refl = new G4OpticalSurface("surface_mppc_refl");
    surface_mppc_refl->SetType(dielectric_dielectric);
    surface_mppc_refl->SetFinish(polished);
    surface_mppc_refl->SetModel(unified);

    // A bare dielectric_dielectric polished surface uses the RINDEX of the two materials
    // (quartz and epoxy) to calculate Fresnel reflection and transmission.
    if (m_kvc_pv) {
      for (const auto mppc_pv : m_mppc_pvs) {
        new G4LogicalBorderSurface("QuartzToMppcRefl", m_kvc_pv, mppc_pv, surface_mppc_refl);
        new G4LogicalBorderSurface("MppcToQuartzRefl", mppc_pv, m_kvc_pv, surface_mppc_refl);
      }
    }
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::DumpMaterialProperties(G4Material* mat)
{
#if DEBUG
  using CLHEP::eV;

  G4cout << "=== Material: " << mat->GetName() << " ===" << G4endl;

  auto mat_prop_table = mat->GetMaterialPropertiesTable();
  if (!mat_prop_table) {
    G4cout << "No material properties table found." << G4endl;
    return;
  }

  const std::vector<G4String> property_names = {"RINDEX", "ABSLENGTH", "REFLECTIVITY"};

  for (const auto& prop : property_names) {
    if (mat_prop_table->ConstPropertyExists(prop)) {
      G4cout << prop << ": " << mat_prop_table->GetConstProperty(prop) << G4endl;
    }
  }

  for (const auto& prop : property_names) {
    G4MaterialPropertyVector* mpv = mat_prop_table->GetProperty(prop);
    if (mpv) {
      G4cout << prop << ":" << G4endl;
      for (size_t i = 0; i < mpv->GetVectorLength(); ++i) {
        G4cout << "Energy: " << mpv->Energy(i) / eV << " eV, "
               << " Value: " << (*mpv)[i] << G4endl;
      }
    }
  }
#else
  (void)mat;
#endif
}
