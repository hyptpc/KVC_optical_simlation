#include "DetectorConstruction.hh"
#include "MPPCSD.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), m_check_overlaps(true)
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
  name = "Phoshorus";
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

}

//_____________________________________________________________________________
void
DetectorConstruction::AddOpticalProperties()
{
  using CLHEP::eV;
  using CLHEP::m;
  using CLHEP::mm;
  
  // +-----------------+
  // | Quartz Property |
  // +-----------------+
  auto quartz_prop = new G4MaterialPropertiesTable();

  std::vector<G4double> photon_energy{
    1.0022 * eV, 1.0045 * eV, 1.0069 * eV, 1.0092 * eV, 1.0115 * eV, 1.0138 * eV,
    1.0162 * eV, 1.0185 * eV, 1.0209 * eV, 1.0232 * eV, 1.0256 * eV, 1.0279 * eV,
    1.0303 * eV, 1.0327 * eV, 1.0351 * eV, 1.0374 * eV, 1.0398 * eV, 1.0422 * eV,
    1.0446 * eV, 1.0470 * eV, 1.0495 * eV, 1.0519 * eV, 1.0543 * eV, 1.0567 * eV,
    1.0592 * eV, 1.0616 * eV, 1.0641 * eV, 1.0665 * eV, 1.0690 * eV, 1.0714 * eV,
    1.0739 * eV, 1.0764 * eV, 1.0789 * eV, 1.0813 * eV, 1.0838 * eV, 1.0863 * eV,
    1.0888 * eV, 1.0914 * eV, 1.0939 * eV, 1.0964 * eV, 1.0989 * eV, 1.1015 * eV,
    1.1040 * eV, 1.1065 * eV, 1.1091 * eV, 1.1116 * eV, 1.1142 * eV, 1.1168 * eV,
    1.1194 * eV, 1.1219 * eV, 1.1245 * eV, 1.1271 * eV, 1.1297 * eV, 1.1323 * eV,
    1.1349 * eV, 1.1375 * eV, 1.1402 * eV, 1.1428 * eV, 1.1454 * eV, 1.1481 * eV,
    1.1507 * eV, 1.1534 * eV, 1.1560 * eV, 1.1587 * eV, 1.1614 * eV, 1.1640 * eV,
    1.1667 * eV, 1.1694 * eV, 1.1721 * eV, 1.1748 * eV, 1.1775 * eV, 1.1802 * eV,
    1.1830 * eV, 1.1857 * eV, 1.1884 * eV, 1.1911 * eV, 1.1939 * eV, 1.1966 * eV,
    1.1994 * eV, 1.2022 * eV, 1.2049 * eV, 1.2077 * eV, 1.2105 * eV, 1.2133 * eV,
    1.2161 * eV, 1.2189 * eV, 1.2217 * eV, 1.2245 * eV, 1.2273 * eV, 1.2302 * eV,
    1.2330 * eV, 1.2359 * eV, 1.2387 * eV, 1.2416 * eV, 1.2444 * eV, 1.2473 * eV,
    1.2502 * eV, 1.2530 * eV, 1.2559 * eV, 1.2588 * eV, 1.2617 * eV, 1.2646 * eV,
    1.2676 * eV, 1.2705 * eV, 1.2734 * eV, 1.2763 * eV, 1.2793 * eV, 1.2822 * eV,
    1.2852 * eV, 1.2882 * eV, 1.2911 * eV, 1.2941 * eV, 1.2971 * eV, 1.3001 * eV,
    1.3031 * eV, 1.3061 * eV, 1.3091 * eV, 1.3121 * eV, 1.3151 * eV, 1.3182 * eV,
    1.3212 * eV, 1.3242 * eV, 1.3273 * eV, 1.3304 * eV, 1.3334 * eV, 1.3365 * eV,
    1.3396 * eV, 1.3427 * eV, 1.3458 * eV, 1.3489 * eV, 1.3520 * eV, 1.3551 * eV,
    1.3582 * eV, 1.3613 * eV, 1.3645 * eV, 1.3676 * eV, 1.3708 * eV, 1.3739 * eV,
    1.3771 * eV, 1.3803 * eV, 1.3835 * eV, 1.3866 * eV, 1.3898 * eV, 1.3930 * eV,
    1.3963 * eV, 1.3995 * eV, 1.4027 * eV, 1.4059 * eV, 1.4092 * eV, 1.4124 * eV,
    1.4157 * eV, 1.4189 * eV, 1.4222 * eV, 1.4255 * eV, 1.4288 * eV, 1.4321 * eV,
    1.4354 * eV, 1.4387 * eV, 1.4420 * eV, 1.4453 * eV, 1.4487 * eV, 1.4520 * eV,
    1.4553 * eV, 1.4587 * eV, 1.4621 * eV, 1.4654 * eV, 1.4688 * eV, 1.4722 * eV,
    1.4756 * eV, 1.4790 * eV, 1.4824 * eV, 1.4858 * eV, 1.4892 * eV, 1.4927 * eV,
    1.4961 * eV, 1.4996 * eV, 1.5030 * eV, 1.5065 * eV, 1.5100 * eV, 1.5134 * eV,
    1.5169 * eV, 1.5204 * eV, 1.5239 * eV, 1.5274 * eV, 1.5310 * eV, 1.5345 * eV,
    1.5380 * eV, 1.5416 * eV, 1.5451 * eV, 1.5487 * eV, 1.5523 * eV, 1.5558 * eV,
    1.5594 * eV, 1.5630 * eV, 1.5666 * eV, 1.5702 * eV, 1.5739 * eV, 1.5775 * eV,
    1.5811 * eV, 1.5848 * eV, 1.5884 * eV, 1.5921 * eV, 1.5958 * eV, 1.5994 * eV,
    1.6031 * eV, 1.6068 * eV, 1.6105 * eV, 1.6142 * eV, 1.6180 * eV, 1.6217 * eV,
    1.6254 * eV, 1.6292 * eV, 1.6329 * eV, 1.6367 * eV, 1.6405 * eV, 1.6442 * eV,
    1.6480 * eV, 1.6518 * eV, 1.6556 * eV, 1.6595 * eV, 1.6633 * eV, 1.6671 * eV,
    1.6710 * eV, 1.6748 * eV, 1.6787 * eV, 1.6825 * eV, 1.6864 * eV, 1.6903 * eV,
    1.6942 * eV, 1.6981 * eV, 1.7020 * eV, 1.7060 * eV, 1.7099 * eV, 1.7138 * eV,
    1.7178 * eV, 1.7217 * eV, 1.7257 * eV, 1.7297 * eV, 1.7337 * eV, 1.7377 * eV,
    1.7417 * eV, 1.7457 * eV, 1.7497 * eV, 1.7537 * eV, 1.7578 * eV, 1.7618 * eV,
    1.7659 * eV, 1.7700 * eV, 1.7741 * eV, 1.7781 * eV, 1.7822 * eV, 1.7863 * eV,
    1.7905 * eV, 1.7946 * eV, 1.7987 * eV, 1.8029 * eV, 1.8070 * eV, 1.8112 * eV,
    1.8154 * eV, 1.8196 * eV, 1.8238 * eV, 1.8280 * eV, 1.8322 * eV, 1.8364 * eV,
    1.8406 * eV, 1.8449 * eV, 1.8491 * eV, 1.8534 * eV, 1.8577 * eV, 1.8619 * eV,
    1.8662 * eV, 1.8705 * eV, 1.8748 * eV, 1.8792 * eV, 1.8835 * eV, 1.8878 * eV,
    1.8922 * eV, 1.8966 * eV, 1.9009 * eV, 1.9053 * eV, 1.9097 * eV, 1.9141 * eV,
    1.9185 * eV, 1.9229 * eV, 1.9274 * eV, 1.9318 * eV, 1.9363 * eV, 1.9407 * eV,
    1.9452 * eV, 1.9497 * eV, 1.9542 * eV, 1.9587 * eV, 1.9632 * eV, 1.9677 * eV,
    1.9723 * eV, 1.9768 * eV, 1.9814 * eV, 1.9859 * eV, 1.9905 * eV, 1.9951 * eV,
    1.9997 * eV, 2.0043 * eV, 2.0089 * eV, 2.0136 * eV, 2.0182 * eV, 2.0229 * eV,
    2.0275 * eV, 2.0322 * eV, 2.0369 * eV, 2.0416 * eV, 2.0463 * eV, 2.0510 * eV,
    2.0557 * eV, 2.0605 * eV, 2.0652 * eV, 2.0700 * eV, 2.0748 * eV, 2.0795 * eV,
    2.0843 * eV, 2.0891 * eV, 2.0939 * eV, 2.0988 * eV, 2.1036 * eV, 2.1085 * eV,
    2.1133 * eV, 2.1182 * eV, 2.1231 * eV, 2.1280 * eV, 2.1329 * eV, 2.1378 * eV,
    2.1427 * eV, 2.1477 * eV, 2.1526 * eV, 2.1576 * eV, 2.1626 * eV, 2.1675 * eV,
    2.1725 * eV, 2.1775 * eV, 2.1826 * eV, 2.1876 * eV, 2.1926 * eV, 2.1977 * eV,
    2.2028 * eV, 2.2078 * eV, 2.2129 * eV, 2.2180 * eV, 2.2231 * eV, 2.2283 * eV,
    2.2334 * eV, 2.2385 * eV, 2.2437 * eV, 2.2489 * eV, 2.2541 * eV, 2.2593 * eV,
    2.2645 * eV, 2.2697 * eV, 2.2749 * eV, 2.2802 * eV, 2.2854 * eV, 2.2907 * eV,
    2.2960 * eV, 2.3013 * eV, 2.3066 * eV, 2.3119 * eV, 2.3172 * eV, 2.3226 * eV,
    2.3279 * eV, 2.3333 * eV, 2.3387 * eV, 2.3440 * eV, 2.3494 * eV, 2.3549 * eV,
    2.3603 * eV, 2.3657 * eV, 2.3712 * eV, 2.3767 * eV, 2.3821 * eV, 2.3876 * eV,
    2.3931 * eV, 2.3986 * eV, 2.4042 * eV, 2.4097 * eV, 2.4153 * eV, 2.4208 * eV,
    2.4264 * eV, 2.4320 * eV, 2.4376 * eV, 2.4432 * eV, 2.4489 * eV, 2.4545 * eV,
    2.4602 * eV, 2.4659 * eV, 2.4715 * eV, 2.4772 * eV, 2.4829 * eV, 2.4887 * eV,
    2.4944 * eV, 2.5002 * eV, 2.5059 * eV, 2.5117 * eV, 2.5175 * eV, 2.5233 * eV,
    2.5291 * eV, 2.5349 * eV, 2.5408 * eV, 2.5466 * eV, 2.5525 * eV, 2.5584 * eV,
    2.5643 * eV, 2.5702 * eV, 2.5761 * eV, 2.5821 * eV, 2.5880 * eV, 2.5940 * eV,
    2.6000 * eV, 2.6060 * eV, 2.6120 * eV, 2.6180 * eV, 2.6240 * eV, 2.6301 * eV,
    2.6361 * eV, 2.6422 * eV, 2.6483 * eV, 2.6544 * eV, 2.6605 * eV, 2.6667 * eV,
    2.6728 * eV, 2.6790 * eV, 2.6851 * eV, 2.6913 * eV, 2.6975 * eV, 2.7037 * eV,
    2.7100 * eV, 2.7162 * eV, 2.7225 * eV, 2.7288 * eV, 2.7351 * eV, 2.7414 * eV,
    2.7477 * eV, 2.7540 * eV, 2.7604 * eV, 2.7667 * eV, 2.7731 * eV, 2.7795 * eV,
    2.7859 * eV, 2.7923 * eV, 2.7988 * eV, 2.8052 * eV, 2.8117 * eV, 2.8182 * eV,
    2.8247 * eV, 2.8312 * eV, 2.8377 * eV, 2.8442 * eV, 2.8508 * eV, 2.8574 * eV,
    2.8640 * eV, 2.8706 * eV, 2.8772 * eV, 2.8838 * eV, 2.8905 * eV, 2.8971 * eV,
    2.9038 * eV, 2.9105 * eV, 2.9172 * eV, 2.9239 * eV, 2.9307 * eV, 2.9374 * eV,
    2.9442 * eV, 2.9510 * eV, 2.9578 * eV, 2.9646 * eV, 2.9714 * eV, 2.9783 * eV,
    2.9852 * eV, 2.9920 * eV, 2.9989 * eV, 3.0058 * eV, 3.0128 * eV, 3.0197 * eV,
    3.0267 * eV, 3.0337 * eV, 3.0406 * eV, 3.0477 * eV, 3.0547 * eV, 3.0617 * eV,
    3.0688 * eV, 3.0759 * eV, 3.0829 * eV, 3.0901 * eV, 3.0972 * eV, 3.1043 * eV,
    3.1115 * eV, 3.1187 * eV, 3.1258 * eV, 3.1330 * eV, 3.1403 * eV, 3.1475 * eV,
    3.1548 * eV, 3.1620 * eV, 3.1693 * eV, 3.1766 * eV, 3.1839 * eV, 3.1913 * eV,
    3.1987 * eV, 3.2060 * eV, 3.2134 * eV, 3.2208 * eV, 3.2282 * eV, 3.2357 * eV,
    3.2431 * eV, 3.2506 * eV, 3.2581 * eV, 3.2656 * eV, 3.2732 * eV, 3.2807 * eV,
    3.2883 * eV, 3.2958 * eV, 3.3034 * eV, 3.3111 * eV, 3.3187 * eV, 3.3263 * eV,
    3.3340 * eV, 3.3417 * eV, 3.3494 * eV, 3.3571 * eV, 3.3649 * eV, 3.3726 * eV,
    3.3804 * eV, 3.3882 * eV, 3.3960 * eV, 3.4038 * eV, 3.4117 * eV, 3.4195 * eV,
    3.4274 * eV, 3.4353 * eV, 3.4432 * eV, 3.4512 * eV, 3.4591 * eV, 3.4671 * eV,
    3.4751 * eV, 3.4831 * eV, 3.4911 * eV, 3.4992 * eV, 3.5072 * eV, 3.5153 * eV,
    3.5234 * eV, 3.5316 * eV, 3.5397 * eV, 3.5479 * eV, 3.5560 * eV, 3.5642 * eV,
    3.5725 * eV, 3.5807 * eV, 3.5889 * eV, 3.5972 * eV, 3.6055 * eV, 3.6138 * eV,
    3.6222 * eV, 3.6305 * eV, 3.6389 * eV, 3.6473 * eV, 3.6557 * eV, 3.6641 * eV,
    3.6725 * eV, 3.6810 * eV, 3.6895 * eV, 3.6980 * eV, 3.7065 * eV, 3.7151 * eV,
    3.7236 * eV, 3.7322 * eV, 3.7408 * eV, 3.7494 * eV, 3.7581 * eV, 3.7667 * eV,
    3.7754 * eV, 3.7841 * eV, 3.7929 * eV, 3.8016 * eV, 3.8104 * eV, 3.8192 * eV,
    3.8279 * eV, 3.8368 * eV, 3.8456 * eV, 3.8545 * eV, 3.8634 * eV, 3.8723 * eV,
    3.8812 * eV, 3.8902 * eV, 3.8991 * eV, 3.9081 * eV, 3.9171 * eV, 3.9261 * eV,
    3.9352 * eV, 3.9443 * eV, 3.9534 * eV, 3.9625 * eV, 3.9716 * eV, 3.9808 * eV,
    3.9899 * eV, 3.9991 * eV, 4.0084 * eV, 4.0176 * eV, 4.0269 * eV, 4.0361 * eV,
    4.0455 * eV, 4.0548 * eV, 4.0641 * eV, 4.0735 * eV, 4.0829 * eV, 4.0923 * eV,
    4.1017 * eV, 4.1112 * eV, 4.1207 * eV, 4.1301 * eV, 4.1397 * eV, 4.1492 * eV,
    4.1588 * eV, 4.1684 * eV, 4.1780 * eV, 4.1876 * eV, 4.1973 * eV, 4.2069 * eV,
    4.2166 * eV, 4.2264 * eV, 4.2361 * eV, 4.2459 * eV, 4.2557 * eV, 4.2655 * eV,
    4.2753 * eV, 4.2852 * eV, 4.2950 * eV, 4.3049 * eV, 4.3149 * eV, 4.3248 * eV,
    4.3348 * eV, 4.3448 * eV, 4.3548 * eV, 4.3648 * eV, 4.3749 * eV, 4.3850 * eV,
    4.3951 * eV, 4.4052 * eV, 4.4154 * eV, 4.4255 * eV, 4.4357 * eV, 4.4460 * eV,
    4.4562 * eV, 4.4665 * eV, 4.4768 * eV, 4.4871 * eV, 4.4974 * eV, 4.5078 * eV,
    4.5182 * eV, 4.5286 * eV, 4.5391 * eV, 4.5495 * eV, 4.5600 * eV, 4.5705 * eV,
    4.5811 * eV, 4.5916 * eV, 4.6022 * eV, 4.6128 * eV, 4.6234 * eV, 4.6341 * eV,
    4.6448 * eV, 4.6555 * eV, 4.6662 * eV, 4.6770 * eV, 4.6878 * eV, 4.6986 * eV,
    4.7094 * eV, 4.7203 * eV, 4.7312 * eV, 4.7420 * eV, 4.7530 * eV, 4.7640 * eV,
    4.7749 * eV, 4.7859 * eV, 4.7970 * eV, 4.8080 * eV, 4.8191 * eV, 4.8302 * eV,
    4.8414 * eV, 4.8525 * eV, 4.8637 * eV, 4.8749 * eV, 4.8862 * eV, 4.8974 * eV,
    4.9087 * eV, 4.9200 * eV, 4.9314 * eV, 4.9427 * eV, 4.9541 * eV, 4.9655 * eV,
    4.9770 * eV, 4.9885 * eV, 5.0000 * eV, 5.0115 * eV, 5.0230 * eV, 5.0346 * eV,
    5.0462 * eV, 5.0579 * eV, 5.0695 * eV, 5.0812 * eV, 5.0929 * eV, 5.1046 * eV,
    5.1164 * eV, 5.1282 * eV, 5.1400 * eV, 5.1519 * eV, 5.1638 * eV, 5.1757 * eV,
    5.1876 * eV, 5.1996 * eV, 5.2115 * eV, 5.2236 * eV, 5.2356 * eV, 5.2477 * eV,
    5.2598 * eV, 5.2719 * eV, 5.2840 * eV, 5.2962 * eV, 5.3084 * eV, 5.3207 * eV,
    5.3329 * eV, 5.3452 * eV, 5.3576 * eV, 5.3699 * eV, 5.3823 * eV, 5.3947 * eV,
    5.4071 * eV, 5.4196 * eV, 5.4321 * eV, 5.4446 * eV, 5.4571 * eV, 5.4697 * eV,
    5.4823 * eV, 5.4950 * eV, 5.5076 * eV, 5.5203 * eV, 5.5331 * eV, 5.5458 * eV,
    5.5586 * eV, 5.5714 * eV, 5.5843 * eV, 5.5972 * eV, 5.6100 * eV, 5.6230 * eV,
    5.6360 * eV, 5.6489 * eV, 5.6619 * eV, 5.6750 * eV, 5.6881 * eV, 5.7012 * eV,
    5.7143 * eV, 5.7275 * eV, 5.7407 * eV, 5.7540 * eV, 5.7672 * eV, 5.7805 * eV,
    5.7938 * eV, 5.8072 * eV, 5.8206 * eV, 5.8340 * eV, 5.8475 * eV, 5.8609 * eV,
    5.8744 * eV, 5.8880 * eV, 5.9016 * eV, 5.9152 * eV, 5.9288 * eV, 5.9425 * eV,
    5.9562 * eV, 5.9699 * eV, 5.9836 * eV, 5.9975 * eV, 6.0113 * eV, 6.0251 * eV,
    6.0390 * eV, 6.0529 * eV, 6.0669 * eV, 6.0809 * eV, 6.0949 * eV, 6.1090 * eV,
    6.1230 * eV, 6.1371 * eV, 6.1513 * eV, 6.1655 * eV, 6.1797 * eV, 6.1939 * eV,
    6.2082 * eV, 6.2225 * eV, 6.2369 * eV, 6.2513 * eV, 6.2657 * eV, 6.2801 * eV,
    6.2946 * eV, 6.3091 * eV, 6.3236 * eV, 6.3382 * eV, 6.3528 * eV, 6.3675 * eV,
    6.3822 * eV, 6.3968 * eV, 6.4116 * eV, 6.4264 * eV, 6.4412 * eV, 6.4560 * eV,
    6.4709 * eV, 6.4859 * eV, 6.5008 * eV, 6.5158 * eV, 6.5308 * eV, 6.5458 * eV,
    6.5609 * eV, 6.5761 * eV, 6.5912 * eV, 6.6064 * eV, 6.6216 * eV, 6.6369 * eV,
    6.6522 * eV, 6.6675 * eV, 6.6829 * eV, 6.6983 * eV, 6.7138 * eV, 6.7292 * eV,
    6.7448 * eV, 6.7603 * eV, 6.7759 * eV, 6.7915 * eV, 6.8072 * eV, 6.8229 * eV,
    6.8386 * eV, 6.8543 * eV, 6.8701 * eV, 6.8860 * eV, 6.9019 * eV, 6.9178 * eV,
    6.9337 * eV, 6.9497 * eV, 6.9657 * eV, 6.9818 * eV, 6.9979 * eV
  };

  std::vector<G4double> refractive_index{
    1.4479, 1.4479, 1.4480, 1.4480, 1.4480, 1.4481, 1.4481, 1.4481,
    1.4482, 1.4482, 1.4482, 1.4482, 1.4483, 1.4483, 1.4483, 1.4484,
    1.4484, 1.4484, 1.4485, 1.4485, 1.4485, 1.4486, 1.4486, 1.4486,
    1.4486, 1.4487, 1.4487, 1.4487, 1.4488, 1.4488, 1.4488, 1.4489,
    1.4489, 1.4489, 1.4490, 1.4490, 1.4490, 1.4490, 1.4491, 1.4491,
    1.4491, 1.4492, 1.4492, 1.4492, 1.4493, 1.4493, 1.4493, 1.4493,
    1.4494, 1.4494, 1.4494, 1.4495, 1.4495, 1.4495, 1.4496, 1.4496,
    1.4496, 1.4496, 1.4497, 1.4497, 1.4497, 1.4498, 1.4498, 1.4498,
    1.4498, 1.4499, 1.4499, 1.4499, 1.4500, 1.4500, 1.4500, 1.4500,
    1.4501, 1.4501, 1.4501, 1.4502, 1.4502, 1.4502, 1.4503, 1.4503,
    1.4503, 1.4503, 1.4504, 1.4504, 1.4504, 1.4505, 1.4505, 1.4505,
    1.4505, 1.4506, 1.4506, 1.4506, 1.4507, 1.4507, 1.4507, 1.4507,
    1.4508, 1.4508, 1.4508, 1.4509, 1.4509, 1.4509, 1.4509, 1.4510,
    1.4510, 1.4510, 1.4511, 1.4511, 1.4511, 1.4511, 1.4512, 1.4512,
    1.4512, 1.4513, 1.4513, 1.4513, 1.4513, 1.4514, 1.4514, 1.4514,
    1.4515, 1.4515, 1.4515, 1.4516, 1.4516, 1.4516, 1.4516, 1.4517,
    1.4517, 1.4517, 1.4518, 1.4518, 1.4518, 1.4518, 1.4519, 1.4519,
    1.4519, 1.4520, 1.4520, 1.4520, 1.4520, 1.4521, 1.4521, 1.4521,
    1.4522, 1.4522, 1.4522, 1.4523, 1.4523, 1.4523, 1.4523, 1.4524,
    1.4524, 1.4524, 1.4525, 1.4525, 1.4525, 1.4526, 1.4526, 1.4526,
    1.4526, 1.4527, 1.4527, 1.4527, 1.4528, 1.4528, 1.4528, 1.4529,
    1.4529, 1.4529, 1.4529, 1.4530, 1.4530, 1.4530, 1.4531, 1.4531,
    1.4531, 1.4532, 1.4532, 1.4532, 1.4533, 1.4533, 1.4533, 1.4533,
    1.4534, 1.4534, 1.4534, 1.4535, 1.4535, 1.4535, 1.4536, 1.4536,
    1.4536, 1.4537, 1.4537, 1.4537, 1.4538, 1.4538, 1.4538, 1.4539,
    1.4539, 1.4539, 1.4540, 1.4540, 1.4540, 1.4541, 1.4541, 1.4541,
    1.4541, 1.4542, 1.4542, 1.4542, 1.4543, 1.4543, 1.4543, 1.4544,
    1.4544, 1.4544, 1.4545, 1.4545, 1.4545, 1.4546, 1.4546, 1.4547,
    1.4547, 1.4547, 1.4548, 1.4548, 1.4548, 1.4549, 1.4549, 1.4549,
    1.4550, 1.4550, 1.4550, 1.4551, 1.4551, 1.4551, 1.4552, 1.4552,
    1.4552, 1.4553, 1.4553, 1.4554, 1.4554, 1.4554, 1.4555, 1.4555,
    1.4555, 1.4556, 1.4556, 1.4556, 1.4557, 1.4557, 1.4558, 1.4558,
    1.4558, 1.4559, 1.4559, 1.4559, 1.4560, 1.4560, 1.4561, 1.4561,
    1.4561, 1.4562, 1.4562, 1.4562, 1.4563, 1.4563, 1.4564, 1.4564,
    1.4564, 1.4565, 1.4565, 1.4566, 1.4566, 1.4566, 1.4567, 1.4567,
    1.4568, 1.4568, 1.4568, 1.4569, 1.4569, 1.4570, 1.4570, 1.4570,
    1.4571, 1.4571, 1.4572, 1.4572, 1.4573, 1.4573, 1.4573, 1.4574,
    1.4574, 1.4575, 1.4575, 1.4576, 1.4576, 1.4576, 1.4577, 1.4577,
    1.4578, 1.4578, 1.4579, 1.4579, 1.4579, 1.4580, 1.4580, 1.4581,
    1.4581, 1.4582, 1.4582, 1.4583, 1.4583, 1.4584, 1.4584, 1.4584,
    1.4585, 1.4585, 1.4586, 1.4586, 1.4587, 1.4587, 1.4588, 1.4588,
    1.4589, 1.4589, 1.4590, 1.4590, 1.4591, 1.4591, 1.4592, 1.4592,
    1.4593, 1.4593, 1.4594, 1.4594, 1.4595, 1.4595, 1.4596, 1.4596,
    1.4597, 1.4597, 1.4598, 1.4598, 1.4599, 1.4599, 1.4600, 1.4600,
    1.4601, 1.4601, 1.4602, 1.4602, 1.4603, 1.4603, 1.4604, 1.4605,
    1.4605, 1.4606, 1.4606, 1.4607, 1.4607, 1.4608, 1.4608, 1.4609,
    1.4610, 1.4610, 1.4611, 1.4611, 1.4612, 1.4612, 1.4613, 1.4614,
    1.4614, 1.4615, 1.4615, 1.4616, 1.4617, 1.4617, 1.4618, 1.4618,
    1.4619, 1.4620, 1.4620, 1.4621, 1.4621, 1.4622, 1.4623, 1.4623,
    1.4624, 1.4624, 1.4625, 1.4626, 1.4626, 1.4627, 1.4628, 1.4628,
    1.4629, 1.4630, 1.4630, 1.4631, 1.4632, 1.4632, 1.4633, 1.4634,
    1.4634, 1.4635, 1.4636, 1.4636, 1.4637, 1.4638, 1.4638, 1.4639,
    1.4640, 1.4640, 1.4641, 1.4642, 1.4643, 1.4643, 1.4644, 1.4645,
    1.4645, 1.4646, 1.4647, 1.4648, 1.4648, 1.4649, 1.4650, 1.4651,
    1.4651, 1.4652, 1.4653, 1.4654, 1.4654, 1.4655, 1.4656, 1.4657,
    1.4657, 1.4658, 1.4659, 1.4660, 1.4661, 1.4661, 1.4662, 1.4663,
    1.4664, 1.4665, 1.4665, 1.4666, 1.4667, 1.4668, 1.4669, 1.4670,
    1.4670, 1.4671, 1.4672, 1.4673, 1.4674, 1.4675, 1.4676, 1.4676,
    1.4677, 1.4678, 1.4679, 1.4680, 1.4681, 1.4682, 1.4683, 1.4684,
    1.4684, 1.4685, 1.4686, 1.4687, 1.4688, 1.4689, 1.4690, 1.4691,
    1.4692, 1.4693, 1.4694, 1.4695, 1.4696, 1.4697, 1.4698, 1.4699,
    1.4700, 1.4701, 1.4702, 1.4703, 1.4704, 1.4705, 1.4706, 1.4707,
    1.4708, 1.4709, 1.4710, 1.4711, 1.4712, 1.4713, 1.4714, 1.4715,
    1.4716, 1.4717, 1.4718, 1.4719, 1.4720, 1.4721, 1.4722, 1.4724,
    1.4725, 1.4726, 1.4727, 1.4728, 1.4729, 1.4730, 1.4731, 1.4733,
    1.4734, 1.4735, 1.4736, 1.4737, 1.4738, 1.4740, 1.4741, 1.4742,
    1.4743, 1.4744, 1.4746, 1.4747, 1.4748, 1.4749, 1.4751, 1.4752,
    1.4753, 1.4754, 1.4756, 1.4757, 1.4758, 1.4759, 1.4761, 1.4762,
    1.4763, 1.4765, 1.4766, 1.4767, 1.4769, 1.4770, 1.4771, 1.4773,
    1.4774, 1.4775, 1.4777, 1.4778, 1.4779, 1.4781, 1.4782, 1.4784,
    1.4785, 1.4787, 1.4788, 1.4789, 1.4791, 1.4792, 1.4794, 1.4795,
    1.4797, 1.4798, 1.4800, 1.4801, 1.4803, 1.4804, 1.4806, 1.4807,
    1.4809, 1.4810, 1.4812, 1.4814, 1.4815, 1.4817, 1.4818, 1.4820,
    1.4822, 1.4823, 1.4825, 1.4827, 1.4828, 1.4830, 1.4832, 1.4833,
    1.4835, 1.4837, 1.4838, 1.4840, 1.4842, 1.4843, 1.4845, 1.4847,
    1.4849, 1.4851, 1.4852, 1.4854, 1.4856, 1.4858, 1.4860, 1.4861,
    1.4863, 1.4865, 1.4867, 1.4869, 1.4871, 1.4873, 1.4875, 1.4876,
    1.4878, 1.4880, 1.4882, 1.4884, 1.4886, 1.4888, 1.4890, 1.4892,
    1.4894, 1.4896, 1.4898, 1.4900, 1.4902, 1.4905, 1.4907, 1.4909,
    1.4911, 1.4913, 1.4915, 1.4917, 1.4919, 1.4922, 1.4924, 1.4926,
    1.4928, 1.4930, 1.4933, 1.4935, 1.4937, 1.4940, 1.4942, 1.4944,
    1.4946, 1.4949, 1.4951, 1.4954, 1.4956, 1.4958, 1.4961, 1.4963,
    1.4966, 1.4968, 1.4970, 1.4973, 1.4975, 1.4978, 1.4981, 1.4983,
    1.4986, 1.4988, 1.4991, 1.4993, 1.4996, 1.4999, 1.5001, 1.5004,
    1.5007, 1.5009, 1.5012, 1.5015, 1.5018, 1.5020, 1.5023, 1.5026,
    1.5029, 1.5032, 1.5034, 1.5037, 1.5040, 1.5043, 1.5046, 1.5049,
    1.5052, 1.5055, 1.5058, 1.5061, 1.5064, 1.5067, 1.5070, 1.5073,
    1.5076, 1.5080, 1.5083, 1.5086, 1.5089, 1.5092, 1.5096, 1.5099,
    1.5102, 1.5105, 1.5109, 1.5112, 1.5115, 1.5119, 1.5122, 1.5126,
    1.5129, 1.5133, 1.5136, 1.5140, 1.5143, 1.5147, 1.5150, 1.5154,
    1.5158, 1.5161, 1.5165, 1.5169, 1.5172, 1.5176, 1.5180, 1.5184,
    1.5188, 1.5192, 1.5195, 1.5199, 1.5203, 1.5207, 1.5211, 1.5215,
    1.5219, 1.5224, 1.5228, 1.5232, 1.5236, 1.5240, 1.5244, 1.5249,
    1.5253, 1.5257, 1.5262, 1.5266, 1.5271, 1.5275, 1.5279, 1.5284,
    1.5289, 1.5293, 1.5298, 1.5302, 1.5307, 1.5312, 1.5317, 1.5321,
    1.5326, 1.5331, 1.5336, 1.5341, 1.5346, 1.5351, 1.5356, 1.5361,
    1.5366, 1.5371, 1.5376, 1.5382, 1.5387, 1.5392, 1.5398, 1.5403,
    1.5408, 1.5414, 1.5419, 1.5425, 1.5430, 1.5436, 1.5442, 1.5448,
    1.5453, 1.5459, 1.5465, 1.5471, 1.5477, 1.5483, 1.5489, 1.5495,
    1.5501, 1.5507, 1.5514, 1.5520, 1.5526, 1.5533, 1.5539, 1.5546,
    1.5552, 1.5559, 1.5566, 1.5572, 1.5579, 1.5586, 1.5593, 1.5600,
    1.5607, 1.5614, 1.5621, 1.5628, 1.5636, 1.5643, 1.5650, 1.5658,
    1.5665, 1.5673, 1.5680, 1.5688, 1.5696, 1.5704, 1.5712, 1.5720,
    1.5728, 1.5736, 1.5744, 1.5752, 1.5761, 1.5769, 1.5778, 1.5786,
    1.5795, 1.5804, 1.5813, 1.5822, 1.5831, 1.5840, 1.5849, 1.5858,
    1.5868, 1.5877, 1.5887, 1.5896, 1.5906
  };

  G4int n_entries = photon_energy.size();
  quartz_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);

  // for absorption_length
  photon_energy = {
    0.9999 * eV, 1.1922 * eV, 1.3776 * eV, 1.5694 * eV, 1.7712 * eV, 1.9680 * eV,
    2.1377 * eV, 2.3393 * eV, 2.5303 * eV, 2.6953 * eV, 2.9520 * eV, 3.0996 * eV,
    3.2627 * eV, 3.6466 * eV, 4.1328 * eV, 4.5920 * eV, 5.1660 * eV, 5.6356 * eV,
    5.9040 * eV, 6.1992 * eV, 6.8880 * eV,
  };
  
  std::vector<G4double> absorption_length{
    100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m,
    100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m,
    100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 9.4029 * m,
    2.3965 * m, 0.5321 * m, 0.0158 * m
  };

  n_entries = photon_energy.size();
  quartz_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);

  m_material_map["QuartzKVC"]->SetMaterialPropertiesTable(quartz_prop);
  

  // +---------+
  // | Aerogel |
  // +---------+
  auto aerogel_prop = new G4MaterialPropertiesTable();
  photon_energy    = { 1.3 * eV, 7.0 * eV };
  refractive_index = { 1.1, 1.1};
  n_entries = photon_energy.size();
  aerogel_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);

  // for absorption_length
  photon_energy = {1.3*eV,1.56*eV,1.68*eV,1.84*eV,2.06*eV,2.26*eV,2.54*eV,2.90*eV,3.10*eV,3.28*eV,3.94*eV,4.94*eV,7.0*eV};
  absorption_length =  {500*mm,128*mm,120*mm,97*mm,77*mm,59*mm,41*mm,26*mm,20*mm,17*mm,8*mm,4*mm,1*mm};
  n_entries = photon_energy.size();
  aerogel_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);

  m_material_map["Aerogel"]->SetMaterialPropertiesTable(aerogel_prop);
  

  // +--------------+
  // | Air Property |
  // +--------------+
  auto air_prop = new G4MaterialPropertiesTable();

  photon_energy    = { 1.3 * eV, 7.0 * eV };
  refractive_index = { 1.0, 1.0};
  n_entries = photon_energy.size();
  air_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  m_material_map["Air"]->SetMaterialPropertiesTable(air_prop);

  
  // +----------------------+
  // | Black sheet Property |
  // +----------------------+
  auto blacksheet_prop = new G4MaterialPropertiesTable();

  photon_energy     = { 1.3 * eV, 7.0 * eV};
  refractive_index  = { 1.6, 1.6 };
  absorption_length = { 1.0e-9 * cm, 1.0e-9 * cm };
  n_entries = photon_energy.size();
  
  blacksheet_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  blacksheet_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Blacksheet"]->SetMaterialPropertiesTable(blacksheet_prop);


  // +-----------------+
  // | Teflon Property |
  // +-----------------+
  auto teflon_prop = new G4MaterialPropertiesTable();

  photon_energy     = { 1.3 * eV, 7.0 * eV};
  refractive_index  = { 1.35, 1.35 };
  absorption_length = { 1.0e-9 * cm, 1.0e-9 * cm };
  n_entries = photon_energy.size();
  
  teflon_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  teflon_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Teflon"]->SetMaterialPropertiesTable(teflon_prop);

  
  // +---------------+
  // | MPPC Property |
  // +---------------+
  auto mppc_prop = new G4MaterialPropertiesTable();

  photon_energy     = { 1.3 * eV, 7.0 * eV};
  refractive_index  = { 1.5, 1.5 };
  absorption_length = { 1.0 * cm, 1.0 * cm };
  n_entries = photon_energy.size();
  
  mppc_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  mppc_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["MPPC"]->SetMaterialPropertiesTable(mppc_prop);

  
  // +-------------------------------+
  // | MPPC surface (Epoxi) Property |
  // +-------------------------------+
  auto epoxi_prop = new G4MaterialPropertiesTable();

  photon_energy     = { 1.3 * eV, 7.0 * eV};
  refractive_index  = { 1.5, 1.5 };
  absorption_length = { 1.0 * cm, 1.0 * cm };
  n_entries = photon_energy.size();
  
  epoxi_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  epoxi_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Epoxi"]->SetMaterialPropertiesTable(epoxi_prop);
  
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructKVC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::eV;
  std::vector<G4double> photon_energy;  
  G4int n_entries;
  
  // G4ThreeVector kvc_size(104.0*mm, 120.0*mm, 10.0*mm);
  G4ThreeVector kvc_size(104.0*mm, 120.0*mm, 20.0*mm);

  G4ThreeVector origin_pos(0.0*mm, 0.0*mm, 0.0*mm);

  // +---------------------+
  // | Mother Volume (Air) |
  // +---------------------+
  auto mother_solid = new G4Box("KvcMotherSolid",
                                kvc_size.x()/2.0 + 20.0*mm,
                                kvc_size.y()/2.0 + 20.0*mm,
                                kvc_size.z()/2.0 + 20.0*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "KvcMotherLV");
  new G4PVPlacement(nullptr, origin_pos, mother_lv,
                    "KvcMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // -- air surface -----
  auto surface_air = new G4OpticalSurface("surface_air");
  surface_air->SetModel(unified);
  surface_air->SetType(dielectric_dielectric);
  surface_air->SetFinish(ground);  
  new G4LogicalSkinSurface("AirSurface", mother_lv, surface_air);

  
  // +----------+
  // | Radiator |
  // +----------+
  auto kvc_solid = new G4Box("KvcSolid",
			     kvc_size.x()/2.0,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0);
  auto kvc_lv = new G4LogicalVolume(kvc_solid, m_material_map["QuartzKVC"], "KvcLV");  
  // auto kvc_lv = new G4LogicalVolume(kvc_solid, m_material_map["Aerogel"], "KvcLV");  
  new G4PVPlacement(nullptr, origin_pos, kvc_lv, "KvcPV",
		    mother_lv, false, 0, m_check_overlaps);
  kvc_lv->SetVisAttributes(G4Colour::Yellow());

  // -- Quartz surface -----
  auto surface_quartz = new G4OpticalSurface("surface_quartz");
  surface_quartz->SetModel(unified);
  surface_quartz->SetType(dielectric_dielectric);
  surface_quartz->SetFinish(polished);

  // auto surface_quartz_prop = new G4MaterialPropertiesTable();
  // photon_energy    = {1.3*eV,1.56*eV,1.61*eV,1.74*eV,1.90*eV,2.05*eV,2.22*eV,2.34*eV,5.42*eV,7*eV};
  // std::vector<G4double> quartz_reflec{0.85, 0.91, 0.93, 0.95, 0.97, 0.98, 1.0, 1.0, 1.0, 1.0};
  // n_entries = photon_energy.size();
  // surface_quartz_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &quartz_reflec[0], n_entries);

  // // photon_energy = {1.77 * eV, 2.07 * eV, 2.48 * eV, 3.10 * eV, 4.13 * eV}; 
  // // n_entries = photon_energy.size();
  // // std::vector<G4double> teflon_specularLobe(n_entries, 0.1);
  // // std::vector<G4double> teflon_specularSpike = {0.05, 0.05, 0.05, 0.05, 0.05};
  // // std::vector<G4double> teflon_backScatter = {0.0, 0.0, 0.0, 0.0, 0.0};
  // // surface_teflon_prop->AddProperty("SPECULARLOBECONSTANT", &photon_energy[0], &teflon_specularLobe[0], n_entries);
  // // surface_teflon_prop->AddProperty("SPECULARSPIKECONSTANT", &photon_energy[0], &teflon_specularSpike[0], n_entries);
  // // surface_teflon_prop->AddProperty("BACKSCATTERCONSTANT", &photon_energy[0], &teflon_backScatter[0], n_entries);
  // surface_quartz->SetMaterialPropertiesTable(surface_quartz_prop);

  new G4LogicalSkinSurface("QuartzSurface", kvc_lv, surface_quartz);

  
  // +------+
  // | MPPC |
  // +------+
  G4ThreeVector mppc_size(6.0*mm, 6.0*mm, 1.0*mm);

  auto mppc_solid = new G4Box("MppcSolid",
			     mppc_size.x()/2.0,
			     mppc_size.y()/2.0,
			     mppc_size.z()/2.0);
  // auto mppc_lv = new G4LogicalVolume(mppc_solid, m_material_map["MPPC"], "MppcLV");
  auto mppc_lv = new G4LogicalVolume(mppc_solid, m_material_map["Epoxi"], "MppcLV");
  
  auto rot = new G4RotationMatrix;
  rot->rotateX(90.0*deg);
  G4int n_mppc = 16;
  G4double offset = 0.0 * mm;
  // for(G4int i=0; i<n_mppc; ++i){
  //   G4ThreeVector pos_up(
  // 		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
  // 		      kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset,
  // 		      0.0*mm);
  //   G4ThreeVector pos_low(
  // 		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
  // 		      -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset,
  // 		      0.0*mm);
  //   new G4PVPlacement(rot, pos_up, mppc_lv, "MppcPV",
  //                     mother_lv, false, i, m_check_overlaps);
  //   new G4PVPlacement(rot, pos_low, mppc_lv, "MppcPV",
  //                     mother_lv, false, i+n_mppc, m_check_overlaps);
  // }

  for(G4int i=0; i<n_mppc; ++i){
    G4ThreeVector pos_up1(
		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
		      kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset,
		      5.0*mm);
    G4ThreeVector pos_up2(
		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
		      kvc_size.y()/2.0 + mppc_size.z()/2.0 + offset,
		      -5.0*mm);
    
    G4ThreeVector pos_low1(
		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
		      -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset,
		      5.0*mm);
    G4ThreeVector pos_low2(
		      -(mppc_size.x() + 0.5*mm) * ((n_mppc-1)/2.0 - i),
		      -kvc_size.y()/2.0 - mppc_size.z()/2.0 - offset,
		      -5.0*mm);
    
    new G4PVPlacement(rot, pos_up1, mppc_lv, "MppcPV",
                      mother_lv, false, i, m_check_overlaps);
    new G4PVPlacement(rot, pos_up2, mppc_lv, "MppcPV",
                      mother_lv, false, i+n_mppc, m_check_overlaps);
    
    new G4PVPlacement(rot, pos_low1, mppc_lv, "MppcPV",
                      mother_lv, false, i+2*n_mppc, m_check_overlaps);
    new G4PVPlacement(rot, pos_low2, mppc_lv, "MppcPV",
                      mother_lv, false, i+3*n_mppc, m_check_overlaps);
  }

  
  mppc_lv->SetVisAttributes(G4Colour::Blue());

  // -- MPPC Surface -----
  auto surface_mppc = new G4OpticalSurface("surface_mppc");
  surface_mppc->SetModel(unified);
  surface_mppc->SetType(dielectric_dielectric);
  // surface_mppc->SetType(dielectric_metal);
  surface_mppc->SetFinish(polished);

  auto surface_mppc_prop = new G4MaterialPropertiesTable();
  photon_energy = {
    1.3811 * eV, 1.3829 * eV, 1.3846 * eV, 1.3864 * eV, 1.3882 * eV, 1.3899 * eV,
    1.3917 * eV, 1.3935 * eV, 1.3953 * eV, 1.3971 * eV, 1.3989 * eV, 1.4007 * eV,
    1.4025 * eV, 1.4043 * eV, 1.4061 * eV, 1.4079 * eV, 1.4098 * eV, 1.4116 * eV,
    1.4134 * eV, 1.4153 * eV, 1.4171 * eV, 1.4190 * eV, 1.4208 * eV, 1.4227 * eV,
    1.4245 * eV, 1.4264 * eV, 1.4283 * eV, 1.4302 * eV, 1.4321 * eV, 1.4339 * eV,
    1.4358 * eV, 1.4377 * eV, 1.4396 * eV, 1.4415 * eV, 1.4435 * eV, 1.4454 * eV,
    1.4473 * eV, 1.4492 * eV, 1.4512 * eV, 1.4531 * eV, 1.4551 * eV, 1.4570 * eV,
    1.4590 * eV, 1.4609 * eV, 1.4629 * eV, 1.4649 * eV, 1.4668 * eV, 1.4688 * eV,
    1.4708 * eV, 1.4728 * eV, 1.4748 * eV, 1.4768 * eV, 1.4788 * eV, 1.4808 * eV,
    1.4828 * eV, 1.4849 * eV, 1.4869 * eV, 1.4889 * eV, 1.4910 * eV, 1.4930 * eV,
    1.4951 * eV, 1.4971 * eV, 1.4992 * eV, 1.5013 * eV, 1.5034 * eV, 1.5054 * eV,
    1.5075 * eV, 1.5096 * eV, 1.5117 * eV, 1.5138 * eV, 1.5159 * eV, 1.5181 * eV,
    1.5202 * eV, 1.5223 * eV, 1.5244 * eV, 1.5266 * eV, 1.5287 * eV, 1.5309 * eV,
    1.5330 * eV, 1.5352 * eV, 1.5374 * eV, 1.5396 * eV, 1.5417 * eV, 1.5439 * eV,
    1.5461 * eV, 1.5483 * eV, 1.5528 * eV, 1.5550 * eV, 1.5572 * eV, 1.5594 * eV,
    1.5617 * eV, 1.5639 * eV, 1.5662 * eV, 1.5684 * eV, 1.5707 * eV, 1.5730 * eV,
    1.5753 * eV, 1.5775 * eV, 1.5798 * eV, 1.5821 * eV, 1.5844 * eV, 1.5868 * eV,
    1.5891 * eV, 1.5914 * eV, 1.5937 * eV, 1.5961 * eV, 1.5984 * eV, 1.6008 * eV,
    1.6031 * eV, 1.6055 * eV, 1.6079 * eV, 1.6103 * eV, 1.6127 * eV, 1.6150 * eV,
    1.6175 * eV, 1.6199 * eV, 1.6223 * eV, 1.6247 * eV, 1.6271 * eV, 1.6296 * eV,
    1.6320 * eV, 1.6345 * eV, 1.6369 * eV, 1.6394 * eV, 1.6419 * eV, 1.6444 * eV,
    1.6469 * eV, 1.6494 * eV, 1.6519 * eV, 1.6544 * eV, 1.6569 * eV, 1.6594 * eV,
    1.6620 * eV, 1.6645 * eV, 1.6671 * eV, 1.6696 * eV, 1.6722 * eV, 1.6748 * eV,
    1.6774 * eV, 1.6800 * eV, 1.6826 * eV, 1.6852 * eV, 1.6878 * eV, 1.6904 * eV,
    1.6930 * eV, 1.6957 * eV, 1.6983 * eV, 1.7010 * eV, 1.7036 * eV, 1.7063 * eV,
    1.7090 * eV, 1.7117 * eV, 1.7144 * eV, 1.7171 * eV, 1.7198 * eV, 1.7225 * eV,
    1.7253 * eV, 1.7280 * eV, 1.7308 * eV, 1.7335 * eV, 1.7363 * eV, 1.7391 * eV,
    1.7419 * eV, 1.7447 * eV, 1.7475 * eV, 1.7503 * eV, 1.7531 * eV, 1.7559 * eV,
    1.7588 * eV, 1.7616 * eV, 1.7645 * eV, 1.7674 * eV, 1.7702 * eV, 1.7760 * eV,
    1.7789 * eV, 1.7818 * eV, 1.7848 * eV, 1.7877 * eV, 1.7907 * eV, 1.7936 * eV,
    1.7966 * eV, 1.7995 * eV, 1.8025 * eV, 1.8055 * eV, 1.8085 * eV, 1.8115 * eV,
    1.8146 * eV, 1.8176 * eV, 1.8207 * eV, 1.8237 * eV, 1.8268 * eV, 1.8298 * eV,
    1.8329 * eV, 1.8360 * eV, 1.8391 * eV, 1.8423 * eV, 1.8454 * eV, 1.8485 * eV,
    1.8517 * eV, 1.8548 * eV, 1.8580 * eV, 1.8612 * eV, 1.8644 * eV, 1.8676 * eV,
    1.8708 * eV, 1.8740 * eV, 1.8773 * eV, 1.8805 * eV, 1.8838 * eV, 1.8870 * eV,
    1.8903 * eV, 1.8936 * eV, 1.8969 * eV, 1.9002 * eV, 1.9036 * eV, 1.9069 * eV,
    1.9102 * eV, 1.9136 * eV, 1.9170 * eV, 1.9204 * eV, 1.9238 * eV, 1.9272 * eV,
    1.9306 * eV, 1.9340 * eV, 1.9375 * eV, 1.9409 * eV, 1.9444 * eV, 1.9479 * eV,
    1.9514 * eV, 1.9549 * eV, 1.9584 * eV, 1.9620 * eV, 1.9655 * eV, 1.9691 * eV,
    1.9726 * eV, 1.9762 * eV, 1.9798 * eV, 1.9834 * eV, 1.9871 * eV, 1.9907 * eV,
    1.9944 * eV, 1.9980 * eV, 2.0017 * eV, 2.0054 * eV, 2.0091 * eV, 2.0128 * eV,
    2.0166 * eV, 2.0203 * eV, 2.0241 * eV, 2.0279 * eV, 2.0316 * eV, 2.0354 * eV,
    2.0393 * eV, 2.0431 * eV, 2.0469 * eV, 2.0508 * eV, 2.0547 * eV, 2.0586 * eV,
    2.0625 * eV, 2.0703 * eV, 2.0743 * eV, 2.0783 * eV, 2.0822 * eV, 2.0862 * eV,
    2.0902 * eV, 2.0943 * eV, 2.0983 * eV, 2.1024 * eV, 2.1064 * eV, 2.1105 * eV,
    2.1146 * eV, 2.1188 * eV, 2.1229 * eV, 2.1271 * eV, 2.1312 * eV, 2.1354 * eV,
    2.1396 * eV, 2.1438 * eV, 2.1481 * eV, 2.1523 * eV, 2.1566 * eV, 2.1609 * eV,
    2.1652 * eV, 2.1695 * eV, 2.1739 * eV, 2.1782 * eV, 2.1826 * eV, 2.1870 * eV,
    2.1914 * eV, 2.1958 * eV, 2.2003 * eV, 2.2047 * eV, 2.2092 * eV, 2.2137 * eV,
    2.2182 * eV, 2.2228 * eV, 2.2273 * eV, 2.2319 * eV, 2.2365 * eV, 2.2411 * eV,
    2.2457 * eV, 2.2504 * eV, 2.2550 * eV, 2.2597 * eV, 2.2644 * eV, 2.2692 * eV,
    2.2739 * eV, 2.2787 * eV, 2.2835 * eV, 2.2883 * eV, 2.2931 * eV, 2.2979 * eV,
    2.3028 * eV, 2.3077 * eV, 2.3126 * eV, 2.3175 * eV, 2.3225 * eV, 2.3275 * eV,
    2.3325 * eV, 2.3375 * eV, 2.3425 * eV, 2.3476 * eV, 2.3527 * eV, 2.3578 * eV,
    2.3629 * eV, 2.3680 * eV, 2.3732 * eV, 2.3784 * eV, 2.3836 * eV, 2.3889 * eV,
    2.3941 * eV, 2.3994 * eV, 2.4047 * eV, 2.4100 * eV, 2.4154 * eV, 2.4208 * eV,
    2.4262 * eV, 2.4316 * eV, 2.4371 * eV, 2.4425 * eV, 2.4480 * eV, 2.4536 * eV,
    2.4591 * eV, 2.4647 * eV, 2.4703 * eV, 2.4759 * eV, 2.4872 * eV, 2.4930 * eV,
    2.4987 * eV, 2.5044 * eV, 2.5102 * eV, 2.5160 * eV, 2.5219 * eV, 2.5277 * eV,
    2.5336 * eV, 2.5396 * eV, 2.5455 * eV, 2.5515 * eV, 2.5575 * eV, 2.5635 * eV,
    2.5696 * eV, 2.5757 * eV, 2.5818 * eV, 2.5879 * eV, 2.5941 * eV, 2.6003 * eV,
    2.6065 * eV, 2.6128 * eV, 2.6191 * eV, 2.6254 * eV, 2.6318 * eV, 2.6382 * eV,
    2.6446 * eV, 2.6510 * eV, 2.6575 * eV, 2.6640 * eV, 2.6706 * eV, 2.6772 * eV,
    2.6838 * eV, 2.6904 * eV, 2.6971 * eV, 2.7038 * eV, 2.7105 * eV, 2.7173 * eV,
    2.7241 * eV, 2.7310 * eV, 2.7379 * eV, 2.7448 * eV, 2.7517 * eV, 2.7587 * eV,
    2.7657 * eV, 2.7728 * eV, 2.7799 * eV, 2.7870 * eV, 2.7942 * eV, 2.8014 * eV,
    2.8086 * eV, 2.8159 * eV, 2.8232 * eV, 2.8305 * eV, 2.8379 * eV, 2.8454 * eV,
    2.8528 * eV, 2.8603 * eV, 2.8679 * eV, 2.8755 * eV, 2.8831 * eV, 2.8908 * eV,
    2.8985 * eV, 2.9062 * eV, 2.9140 * eV, 2.9218 * eV, 2.9297 * eV, 2.9376 * eV,
    2.9456 * eV, 2.9536 * eV, 2.9617 * eV, 2.9697 * eV, 2.9779 * eV, 2.9861 * eV,
    2.9943 * eV, 3.0026 * eV, 3.0109 * eV, 3.0192 * eV, 3.0277 * eV, 3.0361 * eV,
    3.0446 * eV, 3.0532 * eV, 3.0618 * eV, 3.0704 * eV, 3.0791 * eV, 3.0879 * eV,
    3.1055 * eV, 3.1144 * eV, 3.1234 * eV, 3.1324 * eV, 3.1414 * eV, 3.1505 * eV,
    3.1597 * eV, 3.1689 * eV, 3.1782 * eV, 3.1875 * eV, 3.1968 * eV, 3.2063 * eV,
    3.2158 * eV, 3.2253 * eV, 3.2349 * eV, 3.2446 * eV, 3.2543 * eV, 3.2640 * eV,
    3.2739 * eV, 3.2838 * eV, 3.2937 * eV, 3.3037 * eV, 3.3138 * eV, 3.3239 * eV,
    3.3341 * eV, 3.3444 * eV, 3.3547 * eV, 3.3651 * eV, 3.3756 * eV, 3.3861 * eV,
    3.3967 * eV, 3.4073 * eV, 3.4180 * eV, 3.4288 * eV, 3.4396 * eV, 3.4506 * eV,
    3.4616 * eV, 3.4726 * eV, 3.4837 * eV, 3.4949 * eV, 3.5062 * eV, 3.5176 * eV,
    3.5290 * eV, 3.5405 * eV, 3.5521 * eV, 3.5637 * eV, 3.5754 * eV, 3.5872 * eV,
    3.5991 * eV, 3.6111 * eV, 3.6231 * eV, 3.6352 * eV, 3.6474 * eV, 3.6597 * eV,
    3.6721 * eV, 3.6845 * eV, 3.6970 * eV, 3.7097 * eV, 3.7224 * eV, 3.7351 * eV,
    3.7480 * eV, 3.7610 * eV, 3.7741 * eV, 3.7872 * eV, 3.8004 * eV, 3.8138 * eV,
    3.8272 * eV, 3.8407 * eV, 3.8544 * eV, 3.8681 * eV,
  };
  
  std::vector<G4double> mppc_eff{
    0.1071, 0.1071, 0.1078, 0.1078, 0.1086, 0.1086, 0.1093, 0.1100,
    0.1100, 0.1107, 0.1114, 0.1114, 0.1121, 0.1128, 0.1128, 0.1135,
    0.1135, 0.1143, 0.1150, 0.1157, 0.1157, 0.1164, 0.1171, 0.1171,
    0.1178, 0.1185, 0.1185, 0.1192, 0.1200, 0.1207, 0.1207, 0.1214,
    0.1221, 0.1221, 0.1228, 0.1235, 0.1242, 0.1242, 0.1249, 0.1257,
    0.1264, 0.1264, 0.1271, 0.1278, 0.1285, 0.1285, 0.1292, 0.1299,
    0.1306, 0.1314, 0.1314, 0.1321, 0.1328, 0.1335, 0.1342, 0.1342,
    0.1349, 0.1356, 0.1363, 0.1371, 0.1371, 0.1378, 0.1385, 0.1392,
    0.1399, 0.1399, 0.1406, 0.1413, 0.1420, 0.1428, 0.1435, 0.1435,
    0.1442, 0.1449, 0.1456, 0.1463, 0.1470, 0.1477, 0.1477, 0.1485,
    0.1492, 0.1499, 0.1506, 0.1513, 0.1520, 0.1520, 0.1534, 0.1542,
    0.1549, 0.1556, 0.1563, 0.1570, 0.1577, 0.1577, 0.1584, 0.1591,
    0.1599, 0.1606, 0.1613, 0.1620, 0.1627, 0.1634, 0.1641, 0.1648,
    0.1656, 0.1670, 0.1677, 0.1684, 0.1691, 0.1691, 0.1698, 0.1705,
    0.1713, 0.1720, 0.1727, 0.1734, 0.1741, 0.1748, 0.1748, 0.1762,
    0.1777, 0.1784, 0.1791, 0.1798, 0.1805, 0.1812, 0.1819, 0.1834,
    0.1841, 0.1848, 0.1855, 0.1862, 0.1869, 0.1884, 0.1891, 0.1898,
    0.1905, 0.1912, 0.1926, 0.1933, 0.1941, 0.1948, 0.1962, 0.1969,
    0.1976, 0.1983, 0.1990, 0.2005, 0.2012, 0.2019, 0.2026, 0.2040,
    0.2048, 0.2055, 0.2062, 0.2076, 0.2083, 0.2090, 0.2105, 0.2112,
    0.2119, 0.2126, 0.2140, 0.2147, 0.2154, 0.2169, 0.2176, 0.2183,
    0.2190, 0.2204, 0.2211, 0.2219, 0.2233, 0.2247, 0.2261, 0.2268,
    0.2276, 0.2290, 0.2297, 0.2304, 0.2318, 0.2325, 0.2333, 0.2347,
    0.2354, 0.2368, 0.2375, 0.2390, 0.2397, 0.2404, 0.2418, 0.2425,
    0.2439, 0.2454, 0.2461, 0.2475, 0.2482, 0.2496, 0.2504, 0.2518,
    0.2525, 0.2539, 0.2546, 0.2561, 0.2568, 0.2582, 0.2589, 0.2603,
    0.2610, 0.2625, 0.2639, 0.2646, 0.2667, 0.2675, 0.2689, 0.2696,
    0.2710, 0.2717, 0.2732, 0.2746, 0.2753, 0.2767, 0.2789, 0.2796,
    0.2810, 0.2824, 0.2838, 0.2853, 0.2874, 0.2888, 0.2910, 0.2931,
    0.2945, 0.2967, 0.2988, 0.3002, 0.3024, 0.3045, 0.3059, 0.3074,
    0.3088, 0.3102, 0.3116, 0.3131, 0.3145, 0.3159, 0.3173, 0.3188,
    0.3202, 0.3216, 0.3223, 0.3238, 0.3252, 0.3266, 0.3273, 0.3287,
    0.3302, 0.3316, 0.3330, 0.3352, 0.3366, 0.3380, 0.3401, 0.3416,
    0.3430, 0.3444, 0.3458, 0.3473, 0.3494, 0.3508, 0.3523, 0.3537,
    0.3551, 0.3565, 0.3580, 0.3601, 0.3615, 0.3629, 0.3644, 0.3665,
    0.3679, 0.3686, 0.3701, 0.3722, 0.3736, 0.3751, 0.3765, 0.3779,
    0.3793, 0.3808, 0.3822, 0.3836, 0.3857, 0.3872, 0.3886, 0.3900,
    0.3914, 0.3929, 0.3943, 0.3957, 0.3971, 0.3993, 0.4007, 0.4021,
    0.4036, 0.4050, 0.4064, 0.4078, 0.4100, 0.4114, 0.4128, 0.4143,
    0.4157, 0.4178, 0.4192, 0.4207, 0.4221, 0.4235, 0.4249, 0.4264,
    0.4278, 0.4292, 0.4306, 0.4321, 0.4335, 0.4349, 0.4363, 0.4378,
    0.4385, 0.4399, 0.4413, 0.4428, 0.4435, 0.4449, 0.4463, 0.4470,
    0.4485, 0.4499, 0.4506, 0.4520, 0.4527, 0.4542, 0.4556, 0.4563,
    0.4570, 0.4577, 0.4591, 0.4599, 0.4606, 0.4613, 0.4620, 0.4620,
    0.4627, 0.4634, 0.4641, 0.4641, 0.4648, 0.4656, 0.4656, 0.4663,
    0.4670, 0.4670, 0.4677, 0.4684, 0.4684, 0.4691, 0.4691, 0.4698,
    0.4698, 0.4705, 0.4705, 0.4713, 0.4713, 0.4720, 0.4720, 0.4727,
    0.4727, 0.4727, 0.4727, 0.4727, 0.4734, 0.4734, 0.4727, 0.4727,
    0.4727, 0.4727, 0.4720, 0.4720, 0.4713, 0.4705, 0.4698, 0.4691,
    0.4691, 0.4684, 0.4677, 0.4670, 0.4656, 0.4648, 0.4641, 0.4634,
    0.4627, 0.4620, 0.4613, 0.4606, 0.4599, 0.4591, 0.4584, 0.4577,
    0.4570, 0.4563, 0.4556, 0.4542, 0.4534, 0.4520, 0.4513, 0.4499,
    0.4485, 0.4477, 0.4463, 0.4449, 0.4435, 0.4413, 0.4399, 0.4385,
    0.4371, 0.4356, 0.4342, 0.4328, 0.4314, 0.4299, 0.4285, 0.4271,
    0.4242, 0.4221, 0.4207, 0.4192, 0.4171, 0.4150, 0.4121, 0.4100,
    0.4071, 0.4043, 0.4014, 0.3986, 0.3957, 0.3922, 0.3893, 0.3857,
    0.3822, 0.3786, 0.3751, 0.3715, 0.3679, 0.3651, 0.3608, 0.3572,
    0.3537, 0.3501, 0.3466, 0.3423, 0.3387, 0.3352, 0.3323, 0.3287,
    0.3252, 0.3216, 0.3181, 0.3145, 0.3109, 0.3067, 0.3031, 0.2988,
    0.2945, 0.2903, 0.2860, 0.2810, 0.2753, 0.2689, 0.2625, 0.2553,
    0.2482, 0.2404, 0.2318, 0.2233, 0.2147, 0.2069, 0.1983, 0.1898,
    0.1819, 0.1741, 0.1663, 0.1584, 0.1513, 0.1442, 0.1371, 0.1306,
    0.1242, 0.1185, 0.1150, 0.1114, 0.1086, 0.1064
  };

  n_entries = photon_energy.size();
  std::vector<G4double> mppc_reflec(n_entries, 0.05);

  surface_mppc_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &mppc_reflec[0], n_entries);
  // surface_mppc_prop->AddProperty("EFFICIENCY",   &photon_energy[0], &mppc_eff[0], n_entries);
  // surface_mppc->SetMaterialPropertiesTable(surface_mppc_prop);

  // G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  // G4VPhysicalVolume* kvc_pv = store->GetVolume("KvcPV", false);
  // if (!kvc_pv) {
  //   G4cerr << "Error: Could not find KVC physical volume!" << G4endl;
  // } else {
  //   for (auto pv : *store) {
  //     if (pv->GetName() == "MppcPV") {
  // 	new G4LogicalBorderSurface("KVCtoMPPC_Surface", kvc_pv, pv, surface_mppc);
  //     }
  //   }
  // }
 
  // -- resister SD -----
  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  mppc_lv->SetSensitiveDetector(mppcSD);


  // +--------+
  // | Teflon |
  // +--------+
  auto teflon_solid_full = new G4Box("TeflonSolidFull",
			     kvc_size.x()/2.0 + 1.0*mm,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0 + 1.0*mm);
  auto teflon_solid_cut  = new G4Box("TeflonSolidCut",
			     kvc_size.x()/2.0 + 0.99*mm,
			     kvc_size.y()/2.0,
			     kvc_size.z()/2.0 + 0.99*mm);
  
  // G4SubtractionSolid* teflon_solid = new G4SubtractionSolid("TeflonSolid",
  // 							   teflon_solid_full,
  // 							   kvc_solid,
  // 							   nullptr,
  // 							   origin_pos);
  G4SubtractionSolid* teflon_solid = new G4SubtractionSolid("TeflonSolid",
							   teflon_solid_full,
							   teflon_solid_cut,
							   nullptr,
							   origin_pos);
  
  auto teflon_lv = new G4LogicalVolume(teflon_solid, m_material_map["Teflon"], "TeflonLV");
  new G4PVPlacement(nullptr, origin_pos, teflon_lv, "TeflonPV",
		    mother_lv, false, 0, m_check_overlaps);
  teflon_lv->SetVisAttributes(G4Colour::White());

  // -- teflon surface -----
  auto surface_teflon = new G4OpticalSurface("surface_teflon");
  surface_teflon->SetModel(DAVIS);
  surface_teflon->SetType(dielectric_LUTDAVIS);
  surface_teflon->SetFinish(RoughTeflon_LUT);
  // surface_teflon->SetModel(unified);
  // surface_teflon->SetType(dielectric_dielectric);
  // surface_teflon->SetFinish(groundteflonair);
  
  auto surface_teflon_prop = new G4MaterialPropertiesTable();
  photon_energy    = {1.3*eV,1.56*eV,1.61*eV,1.74*eV,1.90*eV,2.05*eV,2.22*eV,2.34*eV,5.42*eV,7*eV};
  std::vector<G4double> teflon_reflec{0.85, 0.91, 0.93, 0.95, 0.97, 0.98, 1.0, 1.0, 1.0, 1.0};
  n_entries = photon_energy.size();
  surface_teflon_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &teflon_reflec[0], n_entries);

  photon_energy = {1.77 * eV, 2.07 * eV, 2.48 * eV, 3.10 * eV, 4.13 * eV}; 
  n_entries = photon_energy.size();
  std::vector<G4double> teflon_specularLobe(n_entries, 0.1);
  std::vector<G4double> teflon_specularSpike = {0.05, 0.05, 0.05, 0.05, 0.05};
  std::vector<G4double> teflon_backScatter = {0.0, 0.0, 0.0, 0.0, 0.0};
  surface_teflon_prop->AddProperty("SPECULARLOBECONSTANT", &photon_energy[0], &teflon_specularLobe[0], n_entries);
  surface_teflon_prop->AddProperty("SPECULARSPIKECONSTANT", &photon_energy[0], &teflon_specularSpike[0], n_entries);
  surface_teflon_prop->AddProperty("BACKSCATTERCONSTANT", &photon_energy[0], &teflon_backScatter[0], n_entries);
  surface_teflon->SetMaterialPropertiesTable(surface_teflon_prop);
  new G4LogicalSkinSurface("TeflonSurface", teflon_lv, surface_teflon);

  
  // +------------+
  // | Blacksheet |
  // +------------+
  auto blacksheet_solid_full = new G4Box("BlacksheetSolidFull",
			     kvc_size.x()/2.0 + 2.0*mm,
			     kvc_size.y()/2.0 + 2.0*mm,
			     kvc_size.z()/2.0 + 2.0*mm);
  auto blacksheet_solid_cut = new G4Box("BlacksheetSolidCut",
			     kvc_size.x()/2.0 + 1.0*mm,
			     kvc_size.y()/2.0 + 1.0*mm,
			     kvc_size.z()/2.0 + 1.0*mm);
  G4SubtractionSolid* blacksheet_solid = new G4SubtractionSolid("BlacksheetSolid",
							   blacksheet_solid_full,
							   blacksheet_solid_cut,
							   nullptr,
							   origin_pos);
  auto blacksheet_lv = new G4LogicalVolume(blacksheet_solid, m_material_map["Blacksheet"], "BlacksheetLV");
  new G4PVPlacement(nullptr, origin_pos, blacksheet_lv, "BlacksheetPV",
		    mother_lv, false, 0, m_check_overlaps);
  blacksheet_lv->SetVisAttributes(G4Colour::Black());
  
  // -- black sheet surface -----
  auto surface_bs = new G4OpticalSurface("surface_bs");
  surface_bs->SetModel(unified);
  surface_bs->SetType(dielectric_metal);
  surface_bs->SetFinish(ground);

  auto surface_bs_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3*eV, 7.0*eV};
  std::vector<G4double> bs_reflec{0.0, 0.0};
  n_entries = photon_energy.size();
  surface_bs_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &bs_reflec[0], n_entries);
  surface_bs->SetMaterialPropertiesTable(surface_bs_prop);
  new G4LogicalSkinSurface("BlackSheetSurface", blacksheet_lv, surface_bs);
  
}


void DetectorConstruction::DumpMaterialProperties(G4Material* mat)
{
  G4cout << "=== Material: " << mat->GetName() << " ===" << G4endl;

  auto matPropTable = mat->GetMaterialPropertiesTable();
  if (!matPropTable) {
    G4cout << "No material properties table found." << G4endl;
    return;
  }

  std::vector<G4String> propertyNames = {"RINDEX", "ABSORPTION", "REFLECTIVITY"};

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
	G4cout << "  Energy: " << mpv->Energy(i) / eV << " eV, "
	       << " Value: " << (*mpv)[i] << G4endl;
      }
    }
  }
}
