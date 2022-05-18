//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B3DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>
#include "G4SDParticleFilter.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4SubtractionSolid.hh"
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true),fPhantomRegion(0), fDetectorArcRegion(0),crytalLogicalVolume(0),patientLogicalVolume(0)
{
  DefineMaterials();
  fPhantomRegion = new G4Region("Phantom"); // declare a region for the phantom volume
  fDetectorArcRegion = new G4Region ("DetectorArc");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ 
  delete fPhantomRegion;
  delete fDetectorArcRegion;
  delete crytalLogicalVolume;
  delete patientLogicalVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  G4Element*  O = man->FindOrBuildElement("O" , isotopes); 
  G4Element* Gd = man->FindOrBuildElement("Gd",isotopes);
  G4Element* S = man->FindOrBuildElement("S",isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);  

  G4Material* edgeMat = man->FindOrBuildMaterial("G4_Mo");
  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  G4Material* Gd2O2S = new G4Material("Gd2O2S",7.32*g/cm3,3); // this is the material type of the detector in CT
  // the density is referred in the link: https://en.wikipedia.org/wiki/Gadolinium_oxysulfide
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);  
  
  Gd2O2S->AddElement(Gd,2);
  Gd2O2S->AddElement(O,2);
  Gd2O2S->AddElement(S,1);
  
  G4Element* Ca = man->FindOrBuildElement("Ca",isotopes);
  G4Element* Cl = man->FindOrBuildElement("Cl",isotopes);
  G4Element* K = man->FindOrBuildElement("K",isotopes);
  G4Element* P = man->FindOrBuildElement("P",isotopes);
  G4Element* H = man->FindOrBuildElement("H",isotopes);
  G4Element* C = man->FindOrBuildElement("C",isotopes);
  G4Element* F = man->FindOrBuildElement("F",isotopes);
  
  G4Material* HPDE = new G4Material("HPDE",0.95*g/cm3,2);
  HPDE->AddElement(C,2);
  HPDE->AddElement(H,4);
  
  G4Material* Acrylic = new G4Material("Acrylic",1.18*g/cm3,3);
  Acrylic->AddElement(C,5);
  Acrylic->AddElement(O,2);
  Acrylic->AddElement(H,8);
  
  G4Material* Delrin = new G4Material("Delrin",1.41*g/cm3,3);
  Delrin->AddElement(C,1);
  Delrin->AddElement(H,2);
  Delrin->AddElement(O,1);
  
  G4Material* PFFE = new G4Material("PTFE",2.2*g/cm3,2);
  PFFE->AddElement(C,2);
  PFFE->AddElement(F,4);
  
  
  G4Material* CaCl2 = new G4Material("CaCl2",2.15*g/cm3,2);
  CaCl2->AddElement(Ca,1);
  CaCl2->AddElement(Cl,2);
  
  G4Material* H2O = new G4Material("H2O",1*g/cm3,2);
  H2O->AddElement(H,2);
  H2O->AddElement(O,1);
  
  G4Material* KHPO4 = new G4Material("KHPO4",2.37*g/cm3,4);
  KHPO4->AddElement(K,1);
  KHPO4->AddElement(H,1);
  KHPO4->AddElement(P,1);
  KHPO4->AddElement(O,4);
  
  G4Material* Ethanol = new G4Material("Ethanol",0.788*g/cm3,3); // C2H5OH
  Ethanol->AddElement(C,2);
  Ethanol->AddElement(H,6);
  Ethanol->AddElement(O,1);
  
  G4Material* Propanol = new G4Material("Propanol",0.803*g/cm3,3); // C3H8O
  Propanol->AddElement(C,3);
  Propanol->AddElement(H,8);
  Propanol->AddElement(O,1);
  
  G4Material* Butanol = new G4Material("Butanol",0.807*g/cm3,3);//C4H10O
  Butanol->AddElement(C,4);
  Butanol->AddElement(H,10);
  Butanol->AddElement(O,10);
  
  
  G4double densityCaClPercent7 = 1.0515*g/cm3;
  G4double densityCaClPercent18 = 1.1534*g/cm3;
  G4double densityCaClPercent23 = 1.2018*g/cm3;
  G4double densityK2HPO4Percent10 = 1.0747*g/cm3;
  G4double densityK2HPO4Percent17 = 1.1485*g/cm3;
  G4double densityK2HPO4Percent29 = 1.2734*g/cm3;
  G4double densityK2HPO4Percent45 = 1.4671*g/cm3;
  
  G4Material* CaCl2SolutionPercent7  = new G4Material("CaCl2SolutionPercent7",densityCaClPercent7,2);
  CaCl2SolutionPercent7->AddMaterial(CaCl2,7*perCent);
  CaCl2SolutionPercent7->AddMaterial(H2O,93*perCent);
  
  G4Material* CaCl2SolutionPercent18 = new G4Material("CaCl2SolutionPercent18",densityCaClPercent18,2);
  CaCl2SolutionPercent18->AddMaterial(CaCl2,18*perCent);
  CaCl2SolutionPercent18->AddMaterial(H2O,82*perCent);
  
  G4Material* CaCl2SolutionPercent23 = new G4Material("CaCl2SolutionPercent23",densityCaClPercent23,2);
  CaCl2SolutionPercent23->AddMaterial(CaCl2,23*perCent);
  CaCl2SolutionPercent23->AddMaterial(H2O,77*perCent);
  
  G4Material* K2HPO4SolutionPercent10 = new G4Material("K2HPO4SolutionPercent10",densityK2HPO4Percent10,2);
  K2HPO4SolutionPercent10->AddMaterial(KHPO4,10*perCent);
  K2HPO4SolutionPercent10->AddMaterial(H2O,90*perCent);
  
  G4Material* K2HPO4SolutionPercent17 = new G4Material("K2HPO4SolutionPercent17",densityK2HPO4Percent17,2);
  K2HPO4SolutionPercent17->AddMaterial(KHPO4,17*perCent);
  K2HPO4SolutionPercent17->AddMaterial(H2O,83*perCent);
  
  G4Material* K2HPO4SolutionPercent29 = new G4Material("K2HPO4SolutionPercent29",densityK2HPO4Percent29,2);
  K2HPO4SolutionPercent29->AddMaterial(KHPO4,29*perCent);
  K2HPO4SolutionPercent29->AddMaterial(H2O,71*perCent);
  
  G4Material* K2HPO4SolutionPercent45 = new G4Material("K2HPO4SolutionPercent45",densityK2HPO4Percent45,2);
  K2HPO4SolutionPercent45->AddMaterial(KHPO4,45*perCent);
  K2HPO4SolutionPercent45->AddMaterial(H2O,55*perCent);
  
  
  
  
  
  
  
  

//   
//   G4Material* CaCl2Percent7 = new G4Material("CaCl2Percent7", densityCaClPercent7,4); // name, density, number of components
//   CaCl2Percent7->AddElement(Ca,GetMassFractionCaCl2(7)["Ca"]*perCent);
//   CaCl2Percent7->AddElement(Cl,GetMassFractionCaCl2(7)["Cl"]*perCent);
//   CaCl2Percent7->AddElement(H, GetMassFractionCaCl2(7)["H"]*perCent);
//   CaCl2Percent7->AddElement(O,GetMassFractionCaCl2(7)["O"]*perCent);
//   cout<<GetMassFractionCaCl2(0.07)["Ca"]<<endl;
//   cout<<GetMassFractionCaCl2(0.07)["Cl"]<<endl;
//   cout<<GetMassFractionCaCl2(0.07)["H"]<<endl;
//   cout<<GetMassFractionCaCl2(0.07)["O"]<<endl;
//   exit(0);
  
  
  
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
//   G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
//   
//   // declare crystal as a MultiFunctionalDetector scorer
//   //  
//   G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
//   G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
//   G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
//   cryst->RegisterPrimitive(primitiv1);
//   
//   cout<<this<<endl;
//  cout<<"GET HERE"<<endl;
//   
//   
//     //particle filter
//   G4String fltName, particleName;
// 
//   G4SDParticleFilter* gammafilter= new G4SDParticleFilter(fltName="gammaFilter",  particleName="gamma"); // we only tally gamma, basically photons here
// 
//    // particle energy filter
//   G4SDParticleWithEnergyFilter* energyfilter = new  G4SDParticleWithEnergyFilter ("energyFilter", 0*keV, 140*keV); // here we tentatively set the minimum energy is 5keV, but we can change it later
//   
//   G4PSPassageCellCurrent* scorer1 = new G4PSPassageCellCurrent("CellCurrent", fCurrent_In); // the second argument here is 1, which means
//   // we only tally the incoming particles
//   G4PSFlatSurfaceCurrent* scorer2 = new G4PSFlatSurfaceCurrent("FlatCurrent", fCurrent_In);
// 
//    cout<<"GET HERE 2"<<endl;
// //   energyfilter->add("gamma");  // Accept only gamma.
// //   energyfilter->show();        // Show accepting condition to stdout.
//   scorer1->SetFilter(gammafilter);
//   scorer2->SetFilter(gammafilter);
// //   scorer1->SetFilter(energyfilter);
// //   scorer2->SetFilter(energyfilter);
//   
//   cryst->RegisterPrimitive(scorer1);
//   cryst->RegisterPrimitive(scorer2); // attach these two scores into the defined G4MultiFunctionalDetector
//   
//   
//   SetSensitiveDetector("CrystalLV",cryst,true); // this is an inherited method from G4VUserDetectorConstruction, the first 
// //   crytalLogicalVolume->SetSensitiveDetector(cryst);
//   // argument is the name of the logical volume, the second argument is the sensitive detector defined
//    cout<<"GET HERE 3"<<endl;
//   // declare patient as a MultiFunctionalDetector scorer
//   //  
//   G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
//   
//   G4SDManager::GetSDMpointer()->AddNewDetector(patient);
//   
//   G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
//   
//   patient->RegisterPrimitive(primitiv2);
//   
//   SetSensitiveDetector("PatientLV",patient,true);
// //   patientLogicalVolume->SetSensitiveDetector(patient);
//   cout<<"GET HERE 4"<<endl;
//   
  
   G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
//   SetSensitiveDetector("CrystalLV",cryst);
  crytalLogicalVolume->SetSensitiveDetector(cryst);
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
//   SetSensitiveDetector("PatientLV",patient);
   patientLogicalVolume->SetSensitiveDetector(patient);
}

void B3DetectorConstruction::BuildAndSetDetectorArray(G4int sourceDetectorPairID)
{
  // Gamma detector Parameters
  //

    
      
   bool waterPhamtom = false;
   bool bodyPhantom = false;
   bool fourRodsPhantom =false;
   bool DICOMPhantom = false;
   bool AirScan = true;

  
  G4double ring_R1 = 1183.148*mm;
  G4double cryst_dX = 2.7795*mm;
  G4double cryst_dY = 1.12*mm;
  G4double cryst_dZ = 1.4*mm; // this is the dimension of the detector
  G4double crystalWrap_dY = 1.25*mm;
  G4double crystalWrap_dX = cryst_dX + crystalWrap_dY-cryst_dY;
  G4double crystalWrap_dZ = 1.7;
  G4double internalFoilThickness = 80*um; // the thickness of the internal foil for anti-scatter
  G4double internalFoilHeight = 31.75*mm;
  G4double edgeFoilThickness = 40*um;
  G4double edgeFoilHeight = 36.5*mm;
  G4double detectorBox_dX = crystalWrap_dX; 
  G4double detectorBox_dY = crystalWrap_dY + internalFoilThickness; 
  G4double detectorBox_dZ = crystalWrap_dZ; // this is defined according to the Philips CT configuration
  G4double dPhi = 2*std::atan(0.5*detectorBox_dY/ring_R1);
  G4double ring_R2 = (detectorBox_dZ+ring_R1)/std::cos(0.5*dPhi);

  G4int nb_cryst = 816; // this is the number of the total photon detectors in the CT
  G4int nb_rings = 1; // here we just consider one array of detectors, so the number of array is 1
  G4double arcAngle = dPhi*nb_cryst; // this the angle to the arch for holding all the CT detector
  
  cout<<"ring_R1="<<ring_R1<<endl;
  cout<<"ring_R2="<<ring_R2<<endl;
  //
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();
//   G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR"); // this regular air, supposed to be the right one
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_Galactic"); // take it as vacuum
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Gd2O2S");
//   G4Material* edgeMat = nist->FindOrBuildMaterial("G4_Mo");
  G4Material* edgeMat = nist->FindOrBuildMaterial("G4_Au");
  G4Material* wrap_mat = nist->FindOrBuildMaterial("G4_Al");      

  G4bool checkOrNotCheckOverlapping = false; // not checking overlapping
//   G4bool checkOrNotCheckOverlapping = fCheckOverlaps; // checking overlapping
  
  //     
  // World
  //
  G4double world_sizeXY = 3.5*ring_R2;
  G4double world_sizeZ  = 100*detector_dZ;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name
                                   
   physWorld = new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOrNotCheckOverlapping);       // checking overlaps 

                 
  //
  // ring
  //
  G4Tubs* solidRing =
//     new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, 0., twopi);
    new G4Tubs("Ring", 0, ring_R2, 0.5*cryst_dX, 0., twopi);
      
  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
                    
  //     
  // define crystal
  //
  G4double dX = cryst_dX, dY = cryst_dY;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
    
    crytalLogicalVolume = logicCryst;
    fDetectorArcRegion->AddRootLogicalVolume(logicCryst); // add detector arc volume to this region
    
    // define wrapping material 

    G4double dXw = crystalWrap_dX;
    G4double dYw = crystalWrap_dY;
    G4double dZw = crystalWrap_dZ;
    G4Box* solidBox = new G4Box("box1",dXw/2,dYw/2, dZw/2);
    G4double deltaDistance = (crystalWrap_dZ-cryst_dZ)/2;
    G4VSolid* substract = new G4SubtractionSolid ("Wrap-Crystal",solidBox,solidCryst,0,G4ThreeVector(0,0,deltaDistance));
    G4LogicalVolume* logicWrap = new G4LogicalVolume (substract,wrap_mat,"WrapLV");
    
    
  
    
    // place crystals and wrapping within a ring 
    
    G4double initialAngle = 3.0/2*pi; // this is the inital position of the detector arc
//     G4double deltaAngle = 2*pi/1320; // this is the rotation angle per rotation step
//     initialAngle = initialAngle + deltaAngle*sourceDetectorPairID;
  //
  ofstream outfile;
  outfile.open("detectorCenterPosition.csv"); // for writing the center position of detectors
  
  
  

  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi+initialAngle-arcAngle/2.0+dPhi/2;
//     G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ+ deltaDistance*2)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
    
    outfile<<position.x()<<","<<position.y()<<","<<position.z()<<endl;
                                    
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      checkOrNotCheckOverlapping);       // checking overlaps 
      
    G4ThreeVector positionw = (ring_R1+0.5*crystalWrap_dZ)*uz;
    G4Transform3D transformw = G4Transform3D(rotm,positionw);
    new G4PVPlacement( transformw,
		       logicWrap,
		       "wrap",
		       logicRing,
		       false,
		       icrys,
		       checkOrNotCheckOverlapping
                     );
  }
  
  // define and place the internal foil
  
  G4double dX_inf = crystalWrap_dX;
  G4double dY_inf = internalFoilThickness;
  G4double dZ_inf = internalFoilHeight;
  G4double dZ_inw  = 0.8*mm;
  G4Box* solidUpperInterFoil = new G4Box("upperInterFoil", dX_inf/2, dY_inf/2, dZ_inf/2);
  G4Box* solidLowerInterFoil = new G4Box("lowerInterFoil", dX_inf/2, dY_inf/2, dZ_inw/2);
  
  
  G4LogicalVolume* logicalUpperInterFoil = 
  new G4LogicalVolume(solidUpperInterFoil,          //its solid
		      edgeMat,           //its material
		      "upperInterFoilLV");        //its name
  G4LogicalVolume* logicalLowerInterFoil =   
  new G4LogicalVolume(solidLowerInterFoil,          //its solid
		      wrap_mat,           //its material
		      "upperInterFoilLV");        //its name
  
    for (G4int icrys = 0; icrys < nb_cryst-1 ; icrys++) {
    G4double phi = icrys*dPhi+initialAngle-arcAngle/2.0+dPhi/2+dPhi/2;
//     G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+ crystalWrap_dZ - 0.5* dZ_inw)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,             //rotation,position
                      logicalLowerInterFoil,            //its logical volume
                      "lowerInterFoil",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      checkOrNotCheckOverlapping);       // checking overlaps 
    
    G4ThreeVector positionw = (ring_R1+ crystalWrap_dZ - dZ_inw-0.5*internalFoilHeight)*uz;
    G4Transform3D transformw = G4Transform3D(rotm,positionw);
    new G4PVPlacement(transformw,             //rotation,position
		  logicalUpperInterFoil,            //its logical volume
		  "upperInterFoil",             //its name
		  logicRing,             //its mother  volume
		  false,                 //no boolean operation
		  icrys,                 //copy number
		  checkOrNotCheckOverlapping);       // checking overlaps 
  }
  
  // define and place the left edge foil
  
  

  G4double dX_ef = crystalWrap_dX;
  G4double dY_ef = edgeFoilThickness;
  G4double dZ_ef = edgeFoilHeight;
  G4Box* solidEdgeFoil = new G4Box("edgeFoil", dX_ef/2, dY_ef/2, dZ_ef/2);
  G4LogicalVolume* logicalEdgeFoil = 
  new G4LogicalVolume(solidEdgeFoil,          //its solid
		      edgeMat,           //its material
		      "edgeFoilLV");        //its name
  G4RotationMatrix rotm  = G4RotationMatrix();
  rotm.rotateY(90*deg); 
  G4double phi = initialAngle-arcAngle/2.0;
  rotm.rotateZ(phi);
  G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
  G4ThreeVector position = (ring_R1+crystalWrap_dZ-dZ_inw-internalFoilHeight+0.5*edgeFoilHeight)*uz;
  G4Transform3D transform = G4Transform3D(rotm,position);
  new G4PVPlacement( transform,
		     logicalEdgeFoil,
		     "edgeFoil",
		     logicRing,
		     false,
		     checkOrNotCheckOverlapping
                   );
  
  // place the right edge foil
  G4RotationMatrix rotmr  = G4RotationMatrix();
  rotmr.rotateY(90*deg); 
  G4double phir = initialAngle+arcAngle/2.0;
  rotmr.rotateZ(phir);
  G4ThreeVector uzr = G4ThreeVector(std::cos(phir),  std::sin(phir),0.);     
  G4ThreeVector positionr = (ring_R1+crystalWrap_dZ-dZ_inw-internalFoilHeight+0.5*edgeFoilHeight)*uzr;
  G4Transform3D transformr = G4Transform3D(rotmr,positionr);
  new G4PVPlacement( transformr,
		    logicalEdgeFoil,
		    "edgeFoil",
		    logicRing,
		    false,
		    checkOrNotCheckOverlapping
		  );
  

  
  
                                                      
  //
  // full detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);
      
  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
                                 
  // 
  // place rings within detector 
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      checkOrNotCheckOverlapping);       // checking overlaps 
  }
                       
  //
  // place detector in world
  //    
  
    G4double deltaAngle = 2*pi/1320; // this is the rotation angle per rotation step
    G4double rotationAngle =  deltaAngle*sourceDetectorPairID; // total rotation angle
    
    G4RotationMatrix rotm1  = G4RotationMatrix();
    rotm1.rotateZ(rotationAngle);  
    G4ThreeVector position1;
    G4double lengthOfSourceDetectorLine = 645*mm;
    G4double x_iso = -lengthOfSourceDetectorLine*sin(deltaAngle*sourceDetectorPairID);
    G4double y_iso = -lengthOfSourceDetectorLine*(1-cos(deltaAngle*sourceDetectorPairID));
    position1.setX(x_iso);
    position1.setY(y_iso);
    G4Transform3D transform1 = G4Transform3D(rotm1,position1);
    

//     G4double x_afterRotation = 
                 
//   new G4PVPlacement(0,                       //no rotation
//                     G4ThreeVector(0,0,0),         //at (0,0,0)
//                     logicDetector,           //its logical volume
//                     "Detector",              //its name
//                     logicWorld,              //its mother  volume
//                     false,                   //no boolean operation
//                     0,                       //copy number
//                     checkOrNotCheckOverlapping);         // checking overlaps 
  
    new G4PVPlacement(transform1,
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOrNotCheckOverlapping);         // checking overlaps
    
    
       
  
   //   G4Material* wrap_mat = nist->FindOrBuildMaterial("G4_Al");  
  G4Material* Ethanol_mat = nist->FindOrBuildMaterial("Ethanol");
  G4Material* Propanol_mat = nist->FindOrBuildMaterial("Propanol");
  G4Material* Butanol_mat = nist->FindOrBuildMaterial("Butanol");
  G4Material* Acetone_mat = nist->FindOrBuildMaterial("G4_ACETONE");
  G4Material* CaCl2SolutionPercent7_mat = nist->FindOrBuildMaterial("CaCl2SolutionPercent7");
  G4Material* CaCl2SolutionPercent18_mat = nist->FindOrBuildMaterial("CaCl2SolutionPercent18");
  G4Material* CaCl2SolutionPercent23_mat = nist->FindOrBuildMaterial("CaCl2SolutionPercent23");
  G4Material* K2HPO4SolutionPercent10_mat = nist->FindOrBuildMaterial("K2HPO4SolutionPercent10");
  G4Material* K2HPO4SolutionPercent17_mat = nist->FindOrBuildMaterial("K2HPO4SolutionPercent17");
  G4Material* K2HPO4SolutionPercent29_mat = nist->FindOrBuildMaterial("K2HPO4SolutionPercent29");
  G4Material* K2HPO4SolutionPercent45_mat = nist->FindOrBuildMaterial("K2HPO4SolutionPercent45");
  

  G4Material* HPDE_mat = nist->FindOrBuildMaterial("HPDE");
  G4Material* Acrylic_mat = nist->FindMaterial("Acrylic");
  G4Material* Delrin_mat = nist->FindOrBuildMaterial("Delrin");
  G4Material* PTFE_mat = nist->FindOrBuildMaterial("PTFE");
   

   
    
  //
  // patient
  //
  G4double patient_radius = 10*cm; // here we define a phatom with diamter 20 cm
  G4double patient_dZ = 50*mm;   // here we define the phantom with thickness 10 mm
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_WATER"); // here we use water
    
   G4VSolid* solidPatient;
   
   if (waterPhamtom)
   {
     solidPatient =
     new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi); // this is the cylinder water phantom
   }
   if (bodyPhantom)
   {
     solidPatient = new G4EllipticalTube("Patient",175.5*mm,129.5*mm,50*mm); // this the eye phantom 
   }
//     
//   G4EllipticalTube* solidPatient = new G4EllipticalTube("Patient",175.5*mm,129.5*mm,50*mm); // this the eye phantom 
  
    G4LogicalVolume* logicPatient =                         
    new G4LogicalVolume(solidPatient,        //its solid
                        patient_mat,         //its material
                        "PatientLV");        //its name
    
   if (AirScan)
   {
	solidPatient =
	new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi); // this is the cylinder water phantom
       G4LogicalVolume* logicPatient =                         
       new G4LogicalVolume(solidPatient,        //its solid
                        default_mat,         //its material
                        "PatientLV");        //its name
       patientLogicalVolume = logicPatient;
      
   }
    
    
    
  G4double inserRadius = 16.5*mm;
  G4Tubs* solidInsert = new G4Tubs("Insert",0.,inserRadius, 0.5*patient_dZ, 0., twopi); // Ethanol EthanolInsert, and other inserts
  
  G4LogicalVolume* logicalEthanol = new G4LogicalVolume(solidInsert,Ethanol_mat,"EthanolInsert");
  G4LogicalVolume* logicalPropanol = new G4LogicalVolume(solidInsert,Propanol_mat,"PropanolInsert");
  G4LogicalVolume* logicalButanol = new G4LogicalVolume(solidInsert,Butanol_mat,"ButanolInsert");
  G4LogicalVolume* logicalAcetone = new G4LogicalVolume(solidInsert,Acetone_mat,"AcetoneInsert");
  G4LogicalVolume* logicalCaCl2SolutionPct7 = new G4LogicalVolume(solidInsert,CaCl2SolutionPercent7_mat,"CaCl2SolutionPCt7");
  G4LogicalVolume* logicalCaCl2SolutionPct18 = new G4LogicalVolume(solidInsert,CaCl2SolutionPercent18_mat,"CaCl2SolutionPCt18");
  G4LogicalVolume* logicalCaCl2SolutionPct23 = new G4LogicalVolume(solidInsert,CaCl2SolutionPercent23_mat,"CaCl2SolutionPCt23");
  G4LogicalVolume* logicalK2HPO4SolutioinPct10 = new G4LogicalVolume(solidInsert,K2HPO4SolutionPercent10_mat,"K2HPO4SolutionPct10");
  G4LogicalVolume* logicalK2HPO4SolutioinPct17 = new G4LogicalVolume(solidInsert,K2HPO4SolutionPercent17_mat,"K2HPO4SolutionPct17");
  G4LogicalVolume* logicalK2HPO4SolutioinPct29 = new G4LogicalVolume(solidInsert,K2HPO4SolutionPercent29_mat,"K2HPO4SolutionPct29");
  G4LogicalVolume* logicalK2HPO4SolutioinPct45 = new G4LogicalVolume(solidInsert,K2HPO4SolutionPercent45_mat,"K2HPO4SolutionPct45");

  
  G4double inserRadius1 = 25*mm;
  G4Tubs* solidInsert1 = new G4Tubs("Insert1",0.,inserRadius1, 0.5*100*mm, 0., twopi); // Ethanol EthanolInsert, and other inserts
  
  G4LogicalVolume* logicalHPDE = new G4LogicalVolume(solidInsert1,HPDE_mat,"HPDE");
  G4LogicalVolume* logicalAcrylic = new G4LogicalVolume(solidInsert1,Acrylic_mat,"Acrylic");
  G4LogicalVolume* logicalDelrin = new G4LogicalVolume(solidInsert1,Delrin_mat,"Delrin");
  G4LogicalVolume* logicalPTFE = new G4LogicalVolume(solidInsert1,PTFE_mat,"PTFE");
  
//   


    if (bodyPhantom)
    {
	G4RotationMatrix rotmInsert  = G4RotationMatrix(); 
	G4ThreeVector positionInsert;
	positionInsert.setX(46.6*mm);
	positionInsert.setY(59.9*mm);
	G4Transform3D transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		    logicalEthanol,
		    "EthanolInsert",
		    logicPatient, // its mother volume
		    false,
		    0,
		    checkOrNotCheckOverlapping
		    );

	positionInsert.setX(79.2*mm);
	positionInsert.setY(6.6*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalPropanol,
		      "PropanolInsert",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(-49.6*mm);
	positionInsert.setY(-73.4*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalButanol,
		      "ButanolInsert",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(-11.2*mm);
	positionInsert.setY(-48.7*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalAcetone,
		      "AcetoneInsert",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(63.8*mm);
	positionInsert.setY(-55.9*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalCaCl2SolutionPct7,
		      "CaCl2Percent7",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(7.4*mm);
	positionInsert.setY(36.6*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalCaCl2SolutionPct18,
		      "CaCl2SolutionPCt18",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	
	positionInsert.setX(-68.9*mm);
	positionInsert.setY(42.2*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalCaCl2SolutionPct23,
		      "CaCl2SolutionPCt23",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(-16.7*mm);
	positionInsert.setY(75.2*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalK2HPO4SolutioinPct10,
		      "K2HPO4SolutionPct10",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	
	positionInsert.setX(10.6*mm);
	positionInsert.setY(-88.6*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalK2HPO4SolutioinPct17,
		      "K2HPO4SolutionPct17",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(-84.4*mm);
	positionInsert.setY(-19.9*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalK2HPO4SolutioinPct29,
		      "K2HPO4SolutionPct29",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
	positionInsert.setX(40.9*mm);
	positionInsert.setY(-15.9*mm);
	transformInsert = G4Transform3D(rotmInsert,positionInsert);
	new G4PVPlacement(transformInsert,
		      logicalK2HPO4SolutioinPct45,
		      "K2HPO4SolutionPct45",
		      logicPatient, // its mother volume
		      false,
		      0,
		      checkOrNotCheckOverlapping
		    );
	
      fPhantomRegion->AddRootLogicalVolume(logicalEthanol);
      fPhantomRegion->AddRootLogicalVolume(logicalPropanol);
      fPhantomRegion->AddRootLogicalVolume(logicalButanol);
      fPhantomRegion->AddRootLogicalVolume(logicalAcetone);
      fPhantomRegion->AddRootLogicalVolume(logicalCaCl2SolutionPct7);
      fPhantomRegion->AddRootLogicalVolume(logicalCaCl2SolutionPct18);
      fPhantomRegion->AddRootLogicalVolume(logicalCaCl2SolutionPct23);
      fPhantomRegion->AddRootLogicalVolume(logicalK2HPO4SolutioinPct10);
      fPhantomRegion->AddRootLogicalVolume(logicalK2HPO4SolutioinPct17);
      fPhantomRegion->AddRootLogicalVolume(logicalK2HPO4SolutioinPct29);
      fPhantomRegion->AddRootLogicalVolume(logicalK2HPO4SolutioinPct45);
    
    }
    
  if (fourRodsPhantom)
  {
      G4RotationMatrix rotmInsert  = G4RotationMatrix(); 
      G4ThreeVector positionInsert;
      positionInsert.setX(-123*mm);
      positionInsert.setY(-645*mm);
      G4Transform3D transformInsert = G4Transform3D(rotmInsert,positionInsert);
      new G4PVPlacement(transformInsert,
			logicalHPDE,
			"HPDE",
			logicWorld,
			false,
			checkOrNotCheckOverlapping
			);
      positionInsert.setX(-34*mm);
      positionInsert.setY(-645*mm);
      transformInsert = G4Transform3D(rotmInsert,positionInsert);
      new G4PVPlacement(transformInsert,
			logicalAcrylic,
			"Acrylic",
			logicWorld,
			false,
			checkOrNotCheckOverlapping
	
			);
      positionInsert.setX(63*mm);
      positionInsert.setY(-645*mm);
      transformInsert = G4Transform3D(rotmInsert,positionInsert);
      new G4PVPlacement(transformInsert,
			logicalDelrin,
			"Delrin",
			logicWorld,
			false,
			checkOrNotCheckOverlapping
	
			);
      positionInsert.setX(153*mm);
      positionInsert.setY(-645*mm);
      transformInsert = G4Transform3D(rotmInsert,positionInsert);
      new G4PVPlacement(transformInsert,
			logicalPTFE,
			"PTFE",
			logicWorld,
			false,
			checkOrNotCheckOverlapping
	
			);
      
    fPhantomRegion->AddRootLogicalVolume(logicalHPDE); // this for the four rods simulatioin
    fPhantomRegion->AddRootLogicalVolume(logicalAcrylic); // this for the four rods simulatioin
    fPhantomRegion->AddRootLogicalVolume(logicalDelrin); // this for the four rods simulatioin
    fPhantomRegion->AddRootLogicalVolume(logicalPTFE); // this for the four rods simulatioin

  }

    
  //place patient in world
  
  if (bodyPhantom || waterPhamtom)
  {
                     
    new G4PVPlacement(0,                       //no rotation
		      G4ThreeVector(0,-645*mm,0),         //at (0,0,0)
		      logicPatient,            //its logical volume
		      "Patient",               //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      checkOrNotCheckOverlapping);         // checking overlaps 
    
     patientLogicalVolume = logicPatient;
     fPhantomRegion->AddRootLogicalVolume(logicPatient); // add patient volume to the region
  }
  
  if (DICOMPhantom)
  {
    
    G4LogicalVolume* logicalPatientFromDICOM = BuildPhantomFromDICOM();
    
    new G4PVPlacement(0,                       //no rotation
		      G4ThreeVector(0,-645*mm,0),         //at (0,0,0)
		      logicalPatientFromDICOM,            //its logical volume
		      "PatientFromDICOM",               //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      checkOrNotCheckOverlapping);         // checking overlap
    fPhantomRegion->AddRootLogicalVolume(logicalPatientFromDICOM); // add patient volume to the region
    patientLogicalVolume = logicalPatientFromDICOM;
  }
    
    
                                          
  // Visualization attributes
    
  G4VisAttributes* ringVisAtt = new G4VisAttributes(G4Colour::Blue());// set the color of the ring
  ringVisAtt->SetVisibility (true);
//   
//   G4VisAttributes* detVisAtt = new G4VisAttributes(G4Colour::Red());// set the color of the detector
//   detVisAtt->SetVisibility (false);
  
//   logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
//   logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());   
  G4VisAttributes* crystVisAtt = new G4VisAttributes(G4Colour::Red());// set the color of the detector
  crystVisAtt->SetVisibility (true);
  
  G4VisAttributes* patientVisAtt = new G4VisAttributes(G4Colour::Green());// set the color of the detector
  patientVisAtt->SetVisibility (true);
  
  G4VisAttributes* edgeFoilVisAtt = new G4VisAttributes(G4Colour::Magenta());
  edgeFoilVisAtt->SetVisibility(true);
  
  G4VisAttributes* wrapVisAtt = new G4VisAttributes(G4Colour::Cyan());
  wrapVisAtt->SetVisibility(true);
  //
  logicRing->SetVisAttributes (ringVisAtt);
//   logicDetector->SetVisAttributes (detVisAtt); 
  logicCryst->SetVisAttributes(crystVisAtt);
  logicPatient->SetVisAttributes(patientVisAtt);
  logicalEdgeFoil->SetVisAttributes(edgeFoilVisAtt);
  logicWrap->SetVisAttributes(wrapVisAtt);
  logicalLowerInterFoil->SetVisAttributes(wrapVisAtt);
  logicalUpperInterFoil->SetVisAttributes(edgeFoilVisAtt);
  logicalEthanol->SetVisAttributes(crystVisAtt); // use the color of crystal
  logicalPropanol->SetVisAttributes(crystVisAtt);
  logicalButanol->SetVisAttributes(crystVisAtt);
  logicalAcetone->SetVisAttributes(crystVisAtt);
  logicalCaCl2SolutionPct7->SetVisAttributes(crystVisAtt);
  logicalCaCl2SolutionPct18->SetVisAttributes(crystVisAtt);
  logicalCaCl2SolutionPct23->SetVisAttributes(crystVisAtt);
  logicalK2HPO4SolutioinPct10->SetVisAttributes(crystVisAtt);
  logicalK2HPO4SolutioinPct17->SetVisAttributes(crystVisAtt);
  logicalK2HPO4SolutioinPct29->SetVisAttributes(crystVisAtt);
  logicalK2HPO4SolutioinPct45->SetVisAttributes(crystVisAtt);
  
  logicalHPDE->SetVisAttributes(crystVisAtt);
  logicalAcrylic->SetVisAttributes(wrapVisAtt);
  logicalDelrin->SetVisAttributes(patientVisAtt);
  logicalPTFE->SetVisAttributes(edgeFoilVisAtt);
  
  
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  
    G4Region* fDetectorArcRegion;
  // add logical volumes to the created target
  

  cout<<"FINISH BUILDING GEOMETRY"<<endl;
  
}


void B3DetectorConstruction::ClearSensitveDetectorStuff()
{
  

//   delete cryst;
//   delete primitiv1;
//   cryst->RegisterPrimitive(primitiv1);
// //   SetSensitiveDetector("CrystalLV",cryst);
//   crytalLogicalVolume->SetSensitiveDetector(cryst);
//   
//   // declare patient as a MultiFunctionalDetector scorer
//   //  
//   G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
//   G4SDManager::GetSDMpointer()->AddNewDetector(patient);
//   G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
//   patient->RegisterPrimitive(primitiv2);
// //   SetSensitiveDetector("PatientLV",patient);
//    patientLogicalVolume->SetSensitiveDetector(patient);
}

std::map<string,G4double> B3DetectorConstruction::GetMassFractionCaCl2(G4double massFractionOfSolution)
{
  G4double massTotal = 100;
  G4double massOfCaCl2 = massTotal*massFractionOfSolution;
  G4double massOfCa = massOfCaCl2*40.0/(40.0+35.5*2);
  G4double massOfCl = massOfCaCl2-massOfCa;
  G4double massOfH2O = massTotal - massOfCaCl2;
  G4double massOfH = massOfH2O*2.0/(2.0+16);
  G4double massOfO = massOfH2O-massOfH;
  std::map<string, G4double> elementMassFractionMap;
  elementMassFractionMap["Ca"] = massOfCa;
  elementMassFractionMap["Cl"] = massOfCl;
  elementMassFractionMap["H"] = massOfH;
  elementMassFractionMap["O"] = massOfO;
  return elementMassFractionMap;
}

std::map<string,G4double> B3DetectorConstruction::GetMassFractionK2HPO4(G4double massFractionOfSolution)
{
  G4double massTotal = 100;
  G4double massOfK2HPO4 = massTotal*massFractionOfSolution;
  G4double massOfK = massOfK2HPO4*(39.0*2)/(39.0*2+1+31+16*4);
  G4double massOfP = massOfK2HPO4*(31.0)/(39.0*2+1+31+16*4);
  G4double massOfH1 = massOfK2HPO4*(1.0)/(39.0*2+1+31+16*4);
  G4double massOfO1 = massOfK2HPO4*(16.0*4)/(39.0*2+1+31+16*4);
  G4double massOfH2O = massTotal - massOfK2HPO4;
  G4double massOfH2 = massOfH2O*(2.0)/(2.0+16);
  G4double massOfO2 = massOfH2O - massOfH2;
  G4double massOfH = massOfH1 + massOfH2;
  G4double massOfO = massOfO1 + massOfO2;
  std::map<string, G4double> elementMassFractionMap;
  elementMassFractionMap["K"] = massOfK;
  elementMassFractionMap["H"] = massOfH;
  elementMassFractionMap["P"] = massOfP;
  elementMassFractionMap["O"] = massOfO;
  return elementMassFractionMap;

}


/// Below we try to implement building phantom geometry based on DICOM images
//  the function BuildPhantomFromDICOM is largely amended from the DICOM.cc file in DICOM example

  #include "DicomRegularDetectorConstruction.hh"
  #include "DicomNestedParamDetectorConstruction.hh"
  #include "DicomPartialDetectorConstruction.hh"
  
  #ifdef G4_DCMTK
  #include "DicomFileMgr.hh"
  #else
  #include "DicomHandler.hh"
  #endif
    
  #include "DicomIntersectVolume.hh"
  #include "QGSP_BIC.hh"
  #include "G4tgrMessenger.hh"
  
  

G4LogicalVolume*  B3DetectorConstruction::BuildPhantomFromDICOM()
{
  
    new G4tgrMessenger;
  char* part = getenv( "DICOM_PARTIAL_PARAM" );
  
  G4bool bPartial = FALSE;
  if( part && G4String(part) == "1" ) {
    bPartial = TRUE;
  }
   DicomDetectorConstruction* theGeometry = 0;
  
  #ifdef G4_DCMTK
    DicomFileMgr* theFileMgr = 0;
  #else
    DicomHandler* dcmHandler = 0;
  #endif
    
    if( !bPartial ){
  #ifdef G4_DCMTK // if we have DCMTK, then we do this
      
      theFileMgr = DicomFileMgr::GetInstance();
      theFileMgr->Convert("Data.dat");
      
  #else // else we do this, which is a regular way
      // Treatment of DICOM images before creating the G4runManager
  //     cout<<"we use regular way to build geometry"<<endl;
  //     exit(0);
      dcmHandler = new DicomHandler;
      dcmHandler->CheckFileFormat();
  //     cout<<"Finishing dcmHandler specification"<<endl;
  //     exit(0);
  #endif
      
      // Initialisation of physics, geometry, primary particles ...
      char* nest = getenv( "DICOM_NESTED_PARAM" );
      if( nest && G4String(nest) == "1" ) {
	G4cout<<"Use DicomNestedParamDetectorConstruction"<<G4endl;
	theGeometry = new DicomNestedParamDetectorConstruction();
      } else {
	G4cout<<"Use DicomRegularDetectorConstruction"<<G4endl;
	theGeometry = new DicomRegularDetectorConstruction();
      }
    } else {
      G4cout<<"Use DicomPartialDetectorConstruction"<<G4endl;
      theGeometry = new DicomPartialDetectorConstruction();
    }    
    
    theGeometry->Construct(); // Construct the geometry, inlcuding everything in DICOM example
    
    G4cout<<"Finishing building geometry"<<G4endl;

    G4cout<<"After initialize geometry in runManager"<<G4endl;
    
    G4VPhysicalVolume* thePhysicalVolume =  theGeometry->GetWorldPhysicalVolume();
    G4LogicalVolume* dicomLogicalVolume = theGeometry->GetWorldLogicalVolume();
    
//     cout<<"GET here "<<endl;
//     if (thePhysicalVolume != NULL)
//     {
//       for (int i = 0;i<100;i++)
//       {
// 	cout<<"get a valued pointer of phyiscal volume"<<endl;
//       }
//     }
//     cout<<"the name is "<<dicomLogicalVolume->GetName()<<endl;

//   for (int i = 0;i<100;i++)
//   {
// //     cout<<"HERE"<<endl;
//     cout<<"Print test is "<<theGeometry->PrintTest()<<endl;
// //     cout<<"the name is "<<dicomLogicalVolume->GetName()<<endl;
//     
//   }
    G4LogicalVolume* theLogicalVolume = thePhysicalVolume->GetLogicalVolume();
     
    
     new DicomIntersectVolume();
      
	if( !bPartial ) {
  #ifdef G4_DCMTK
      delete theFileMgr;
  #else
      delete dcmHandler;
  #endif
    }
    
    return dicomLogicalVolume;


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
