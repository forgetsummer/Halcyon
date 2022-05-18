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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "TruebeamHead.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include "G4TwoVector.hh"
#include <fstream>
#include "ReadLinacRTInformation.hh"
#include "ArgumentInterpreter.hh"
#include "G4UserLimits.hh"
using namespace std;


TruebeamHead::TruebeamHead()
{
  DefineMaterials();
}

TruebeamHead::~TruebeamHead()
{
}


void TruebeamHead::DefineMaterials()
{

      G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
      G4Material* Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
      G4Material* Ni = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
      G4Material* tungstenAlloy = new G4Material("tungstenAlloy",18*g/cm3,3);// define the Tungsten Fe Nickel alloy
      G4double fw = 95.0*perCent; // suppose the Tungsten mass fraction is 75%
      G4double fFe = 2.645*perCent; // this is calcualted assuming mass fraction of Cu and Ni are equal
      G4double fNi = 1-fFe-fw; // mass fraction information is from https://www.eaglealloys.com/wp-content/uploads/2017/05/EA-Tungsten-Alloy-Data-Sheet.pdf
      tungstenAlloy->AddMaterial(W,fw);
      tungstenAlloy->AddMaterial(Fe,fFe);
      tungstenAlloy->AddMaterial(Ni,fNi);
      
      
      G4Material* C = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");
      G4Material* Mn = G4NistManager::Instance()->FindOrBuildMaterial("G4_Mn");
      G4Material* S = G4NistManager::Instance()->FindOrBuildMaterial("G4_S");
      G4Material* P = G4NistManager::Instance()->FindOrBuildMaterial("G4_P");
      
      G4double fC = 5*perCent;
      G4double fMn = 0.9*perCent;
      G4double fS = 0.05*perCent;
      G4double fP = 0.04*perCent;
      G4double fFe1 = 1- fC - fMn - fS - fP; // the mass fraction information is from https://www.azom.com/article.aspx?ArticleID=9153
      
      
//       cout<<"The fNi = "<<fNi<<endl;
//       cout<<"The fFe1 = "<<fFe1<<endl;
      
      G4Material* steel = new G4Material("steel",7.87*g/cm3,5);// define the Tungsten Fe Nickel alloy
      steel->AddMaterial(C,fC);
      steel->AddMaterial(Mn, fMn);
      steel->AddMaterial(S, fS);
      steel->AddMaterial(P, fP);
      steel->AddMaterial(Fe, fFe1);
      
      
  
    
}



void TruebeamHead::GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume)
{
  // create the accelerator-world box
  
    PVWorld = physicalVolume;
    

    G4double ssd = 100*cm; // this is the SSD distance, here we take it 100 cm
    SetSSD(ssd);
    SetFieldSize(10*cm);
    
//     Jaw1X();
//     Jaw2X();
//     Jaw1Y();
//     Jaw2Y();
//     BasePlate();
    
    MLC();
    
//     PhaseSpacePlane(733*mm);
//     SourceMLCLeakBlock(500*mm);
    
    G4LogicalVolume* testLV = PVWorld->GetLogicalVolume();
    
}




void TruebeamHead::Jaw1X()
{
    G4String name="Jaw1X";
    G4ThreeVector centerPos;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
    
    G4double topSurfaceXJawToISO = 63.39*cm;
    G4double XJawWith = 5.3*2.54*cm;
    G4double XJawLength = 8.6*2.54*cm;
    G4double XJawHeight = 7.8*cm;
    
    centerPos.set(XJawWith/2.0, 0 ,topSurfaceXJawToISO - XJawHeight/2.0);
    
    G4ThreeVector centerPosForRotate ;
    centerPosForRotate.set(centerPos.getX(),centerPos.getY(), 100*cm - centerPos.getZ());
    
    G4Box* box = new G4Box(name+"Box", XJawWith/2.0, XJawLength/2.0, XJawHeight/2.0);
    G4LogicalVolume* logVol = new G4LogicalVolume(box, tungstenAlloy, name+"LV", 0, 0, 0);
    
    G4double theta = atan(fieldSize/2/SSD);
    
    G4RotationMatrix *cRotation=new G4RotationMatrix();
    cRotation->rotateY(theta);
    
    G4ThreeVector positionAfterRotation = GetCenterPositionOfRotatedJaw(name,centerPosForRotate,theta);
    cout<<"SSD = "<<SSD<<endl;
    cout<<"topSurfaceXJawToISO = "<<topSurfaceXJawToISO<<endl;
    cout<<"XJawHeight/2 = "<<XJawHeight/2<<endl;
    
    cout<<"SSD - topSurfaceXJawToISO + XJawHeight/2.0 = "<<SSD - topSurfaceXJawToISO + XJawHeight/2.0<<endl;
   
    cout<<"theta = "<<theta/deg<<endl;
    cout<<"centerPos = "<<centerPos <<endl;
    positionAfterRotation.setZ(100*cm - positionAfterRotation.getZ());
    cout<<"positionAfterRotation is "<<positionAfterRotation.getX()<<" "<<positionAfterRotation.getY()<<
    " "<<positionAfterRotation.getZ()<<endl;
//     exit(0);
    phVol1X= new G4PVPlacement(cRotation, positionAfterRotation, name+"PV", logVol, PVWorld, false, 0);
    // Region for cuts
    G4Region *regVol;
    regVol= new G4Region(name+"R");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm);
    regVol->SetProductionCuts(cuts);
    logVol->SetRegion(regVol);
    regVol->AddRootLogicalVolume(logVol);

    // Visibility
    G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
    simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
    logVol->SetVisAttributes(simpleAlSVisAtt);
}

void TruebeamHead::Jaw2X()
{
    G4String name="Jaw2X";
    G4ThreeVector centerPos;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
    	  
    G4double topSurfaceXJawToISO = 63.39*cm;
    G4double XJawWith = 5.3*2.54*cm;
    G4double XJawLength = 8.6*2.54*cm;
    G4double XJawHeight = 7.8*cm;
    
    centerPos.set(-1*XJawWith/2.0, 0,topSurfaceXJawToISO - XJawHeight/2.0);
    
    G4ThreeVector centerPosForRotate ;
    centerPosForRotate.set(centerPos.getX(),centerPos.getY(), 100*cm - centerPos.getZ());
    
    G4Box* box = new G4Box(name+"Box", XJawWith/2.0, XJawLength/2.0, XJawHeight/2.0);
    G4LogicalVolume* logVol = new G4LogicalVolume(box, tungstenAlloy, name+"LV", 0, 0, 0);
    
    G4double theta = atan(fieldSize/2/SSD);
    G4RotationMatrix *cRotation=new G4RotationMatrix();
    cRotation->rotateY(-theta);
    G4ThreeVector positionAfterRotation = GetCenterPositionOfRotatedJaw(name,centerPosForRotate,theta);
    
    cout<<"theta = "<<theta/deg<<endl;
    cout<<"centerPos = "<<centerPos <<endl;
    positionAfterRotation.setZ(100*cm - positionAfterRotation.getZ());
    cout<<"positionAfterRotation is "<<positionAfterRotation.getX()<<" "<<positionAfterRotation.getY()<<
    " "<<positionAfterRotation.getZ()<<endl;
//     exit(0);
    phVol1X= new G4PVPlacement(cRotation, positionAfterRotation, name+"PV", logVol, PVWorld, false, 0);

    // Region for cuts
    G4Region *regVol;
    regVol= new G4Region(name+"R");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm);
    regVol->SetProductionCuts(cuts);
    logVol->SetRegion(regVol);
    regVol->AddRootLogicalVolume(logVol);

    // Visibility
    G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
    simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
    logVol->SetVisAttributes(simpleAlSVisAtt);

}

void TruebeamHead::Jaw1Y()
{
    G4String name="Jaw1Y";
    G4ThreeVector centerPos;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
    
    G4double topSurfaceYJawToISO = 72.11*cm;
    G4double YJawWith =7.4*2.54*cm; 
    G4double YJawLength = 4.7*2.54*cm;
    G4double YJawHeight = 7.77*cm;
    
    centerPos.set(0,YJawLength/2.0, topSurfaceYJawToISO - YJawHeight/2.0);
       
    G4ThreeVector centerPosForRotate ;
    centerPosForRotate.set(centerPos.getX(),centerPos.getY(), 100*cm - centerPos.getZ());
    
    G4Box* box = new G4Box(name+"Box", YJawWith/2.0, YJawLength/2.0, YJawHeight/2.0);
    
    G4LogicalVolume* logVol = new G4LogicalVolume(box, tungstenAlloy, name+"LV", 0, 0, 0);
    
    G4double theta = atan(fieldSize/2/SSD);
    G4RotationMatrix *cRotation=new G4RotationMatrix();
    cRotation->rotateX(-theta);
    
    G4ThreeVector positionAfterRotation = GetCenterPositionOfRotatedJaw(name,centerPosForRotate,theta);
    cout<<"positionAfterRotation.getZ() = "<<positionAfterRotation.getZ()<<endl;
    positionAfterRotation.setZ(100*cm - positionAfterRotation.getZ());
//     positionAfterRotation = centerPos;
    
    cout<<"centerPos = "<<centerPos <<endl;
    
    cout<<"positionAfterRotation is "<<positionAfterRotation.getX()<<" "<<positionAfterRotation.getY()<<
    " "<<positionAfterRotation.getZ()<<endl;
//     exit(0);
    phVol1X= new G4PVPlacement(cRotation, positionAfterRotation, name+"PV", logVol, PVWorld, false, 0);

    // Region for cuts
    G4Region *regVol;
    regVol= new G4Region(name+"R");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm);
    regVol->SetProductionCuts(cuts);
    logVol->SetRegion(regVol);
    regVol->AddRootLogicalVolume(logVol);

    // Visibility
    G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
    simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
    logVol->SetVisAttributes(simpleAlSVisAtt);
}

void TruebeamHead::Jaw2Y()
{
    G4String name="Jaw2Y";
    G4ThreeVector centerPos, halfSize;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
    G4double topSurfaceYJawToISO = 72.11*cm;
    G4double YJawWith =7.4*2.54*cm; 
    G4double YJawLength = 4.7*2.54*cm;
    G4double YJawHeight = 7.77*cm;
    
    centerPos.set(0,-YJawLength/2.0, topSurfaceYJawToISO - YJawHeight/2.0);
    
    G4ThreeVector centerPosForRotate ;
    centerPosForRotate.set(centerPos.getX(),centerPos.getY(), 100*cm - centerPos.getZ());
    
    G4Box* box = new G4Box(name+"Box", YJawWith/2.0, YJawLength/2.0, YJawHeight/2.0);
    G4LogicalVolume* logVol = new G4LogicalVolume(box, tungstenAlloy, name+"LV", 0, 0, 0);
    
    G4double theta = atan(fieldSize/2/SSD);
    G4RotationMatrix *cRotation=new G4RotationMatrix();
    cRotation->rotateX(theta);
    G4ThreeVector positionAfterRotation = GetCenterPositionOfRotatedJaw(name,centerPosForRotate,theta);
    positionAfterRotation.setZ(100*cm - positionAfterRotation.getZ());
    cout<<"centerPos = "<<centerPos <<endl;
    
    cout<<"positionAfterRotation is "<<positionAfterRotation.getX()<<" "<<positionAfterRotation.getY()<<
    " "<<positionAfterRotation.getZ()<<endl;
//     exit(0);
    phVol1X= new G4PVPlacement(cRotation, positionAfterRotation, name+"PV", logVol, PVWorld, false, 0);

    // Region for cuts
    G4Region *regVol;
    regVol= new G4Region(name+"R");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm);
    regVol->SetProductionCuts(cuts);
    logVol->SetRegion(regVol);
    regVol->AddRootLogicalVolume(logVol);

    // Visibility
    G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Magenta());
    simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
    logVol->SetVisAttributes(simpleAlSVisAtt);

}






void TruebeamHead::PhaseSpacePlane(G4double zPosition)
{
        // added a plane for phase space file   
  
    //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR"); // this regular air, supposed to be the right one
    
  G4Box* phaseSpacePlane = new G4Box("phaseSpacePlane", 1000*mm/2, 80*um/2, 1000*mm/2);
  G4LogicalVolume* phaseSpacePlane_log = 
  new G4LogicalVolume(phaseSpacePlane,          //its solid
		      default_mat,           //its material
		      "phaseSpacePlane");        //its name
  

  G4RotationMatrix rotm1  = G4RotationMatrix();
  rotm1.rotateX(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis  
  
  
  G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
  ReadLinacRTInformation readHlcy;
  readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information
  
  G4double gantryRtn = readHlcy.GetGantryRtnInformation()*deg ;
  G4double collRtn = readHlcy.GetCollimationRtnInformation()*deg;
//   rotm1.rotateX(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
//   rotm1.rotateZ(-180*deg);// Get back to the orginal rotation from current angle
  
//   rotm1.rotateY(gantryRtn);
//   rotm1.rotateZ(collRtn);
  rotm1.rotateZ(collRtn); // rotate with Z axis first
  rotm1.rotateY(gantryRtn);
  
  
  G4ThreeVector oldPosition = G4ThreeVector(0,0,zPosition);
  G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn, collRtn);
  
//   cout<<"xr = "<<newPosition.getX()<<" yr = "<<newPosition.getY()<<" zr = "<<newPosition.getZ()<<endl;
//   exit(0);
  
  G4Transform3D transform1 = G4Transform3D(rotm1,newPosition);
  new G4PVPlacement(transform1,"phaseSpacePlane",phaseSpacePlane_log,PVWorld,false,0,false);
  
  G4VisAttributes* simpleAlSVisAtt;
  simpleAlSVisAtt= new G4VisAttributes(G4Colour::Gray());
  simpleAlSVisAtt->SetVisibility(true);
  phaseSpacePlane_log->SetVisAttributes(simpleAlSVisAtt);

}


void TruebeamHead::SourceMLCLeakBlock(G4double zPosition)
{
      // add a leaking block for the MLC
    G4double innerRadius = 10*2.54*cm;
    G4double outerRadius = 16*2.54*cm;
    G4double thickNess = 10*cm;
    
    
    G4Box* block1 = new G4Box("block1", 1000*mm/2, thickNess, 300*mm/2);
    
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
    G4LogicalVolume* block1_log = new G4LogicalVolume(block1,tungstenAlloy,"block1_log");
    G4LogicalVolume* block2_log = new G4LogicalVolume(block1,tungstenAlloy,"block2_log");
    
    G4RotationMatrix rotm1  = G4RotationMatrix();
    rotm1.rotateX(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis  
    
    G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
    ReadLinacRTInformation readHlcy;
    readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information
    
    G4double gantryRtn = readHlcy.GetGantryRtnInformation()*deg ;
    G4double collRtn = readHlcy.GetCollimationRtnInformation()*deg;
  //   rotm1.rotateX(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
  //   rotm1.rotateZ(-180*deg);// Get back to the orginal rotation from current angle
    
  //   rotm1.rotateY(gantryRtn);
  //   rotm1.rotateZ(collRtn);
    rotm1.rotateZ(collRtn); // rotate with Z axis first
    rotm1.rotateY(gantryRtn);
    
    
    G4ThreeVector oldPosition = G4ThreeVector(0,250*mm,zPosition);
    G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn, collRtn);
    
  //   cout<<"xr = "<<newPosition.getX()<<" yr = "<<newPosition.getY()<<" zr = "<<newPosition.getZ()<<endl;
  //   exit(0);
    
    G4Transform3D transform1 = G4Transform3D(rotm1,newPosition);
    
    G4ThreeVector oldPosition2 = G4ThreeVector(0,-250*mm,zPosition);
    G4ThreeVector newPosition2 = GetCenterPositionOfRotatedMLC(oldPosition2,gantryRtn, collRtn);
    
    G4Transform3D transform2 = G4Transform3D(rotm1,newPosition2);

    
    new G4PVPlacement(transform1,"block1",block1_log,PVWorld,false,0,false);
    new G4PVPlacement(transform2,"block2",block2_log,PVWorld,false,0,false);
    
    G4VisAttributes*  simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
    simpleAlSVisAtt->SetVisibility(true);
    block1_log->SetVisAttributes(simpleAlSVisAtt);
    block2_log->SetVisAttributes(simpleAlSVisAtt);
    
    
            // Region for cuts
    G4Region *regVol;
    regVol= new G4Region("BlockRegion");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm); // set the production cut of MLC
    regVol->SetProductionCuts(cuts);
    block1_log->SetRegion(regVol);
    block2_log->SetRegion(regVol);
    

}


void TruebeamHead::BasePlate()
{
    G4double innerRadius = 6*2.54*cm;
    G4double outerRadius = 11.9*2.54*cm;
    G4double thickNess = 0.6*2.54*cm;
    G4Tubs* basePlate = new G4Tubs("basePlate",innerRadius,outerRadius, 0.5*thickNess, 0., twopi);
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* steel = nist->FindOrBuildMaterial("steel");
    G4LogicalVolume* basePlate_log = new G4LogicalVolume(basePlate,steel, "basePlate_log");
    G4RotationMatrix rotm = G4RotationMatrix();
    
    G4double topSurfaceBasePlaceToISO  = 53.3*cm;
    G4double z_basePlate = topSurfaceBasePlaceToISO;
    G4Transform3D transform = G4Transform3D(rotm,G4ThreeVector(0,0,z_basePlate));
    
    new G4PVPlacement(transform,"baseplate",basePlate_log,PVWorld,false,0,false);
    
    G4VisAttributes*  simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
    simpleAlSVisAtt->SetVisibility(true);
    basePlate_log->SetVisAttributes(simpleAlSVisAtt);
    
}


G4ThreeVector TruebeamHead::GetCenterPositionOfRotatedJaw(G4String jawName,G4ThreeVector initialPosition, G4double theta)
{
  // this function is getting the center position of the jaw after rotation 
  // here we define Xjaw1 is the jaw toward the positive direction of x axis
  // XJaw2 is the jaw toward the negative direction of axis
  // when xjaw rotates, the rotation angle is with respect to y axis
  // when yjaw rotates, the rotation angle is with respect to x axis
  
  std::vector<std::vector<double> > transitionMatrix;
  transitionMatrix =   MatrixInitialization(3,3);
  if (jawName == "Jaw1X")
  {
    
    transitionMatrix[0][0] = cos(theta);
    transitionMatrix[0][1] = 0;
    transitionMatrix[0][2] = sin(theta);
    transitionMatrix[1][0] = 0;
    transitionMatrix[1][1] = 1;
    transitionMatrix[1][2] = 0;
    transitionMatrix[2][0] = -sin(theta);
    transitionMatrix[2][1] = 0;
    transitionMatrix[2][2] = cos(theta);

  }
  if (jawName == "Jaw2X")
  {
    transitionMatrix[0][0] = cos(theta);
    transitionMatrix[0][1] = 0;
    transitionMatrix[0][2] = -sin(theta);
    transitionMatrix[1][0] = 0;
    transitionMatrix[1][1] = 1;
    transitionMatrix[1][2] = 0;
    transitionMatrix[2][0] = sin(theta);
    transitionMatrix[2][1] = 0;
    transitionMatrix[2][2] = cos(theta);
    
  }
  if (jawName == "Jaw1Y")
  {
    transitionMatrix[0][0] = 1;
    transitionMatrix[0][1] = 0;
    transitionMatrix[0][2] = 0;
    transitionMatrix[1][0] = 0;
    transitionMatrix[1][1] = cos(theta);
    transitionMatrix[1][2] = sin(theta);
    transitionMatrix[2][0] = 0;
    transitionMatrix[2][1] = -sin(theta);
    transitionMatrix[2][2] = cos(theta);
    

    
    
  }
  if (jawName == "Jaw2Y")
  {
    transitionMatrix[0][0] = 1;
    transitionMatrix[0][1] = 0;
    transitionMatrix[0][2] = 0;
    transitionMatrix[1][0] = 0;
    transitionMatrix[1][1] = cos(theta);
    transitionMatrix[1][2] = -sin(theta);
    transitionMatrix[2][0] = 0;
    transitionMatrix[2][1] = sin(theta);
    transitionMatrix[2][2] = cos(theta);
    
    

  }
  
//   PrintMatrix(transitionMatrix);
  
  
  
  std::vector<std::vector<double> > oldPosition;
  oldPosition = MatrixInitialization(3,1);
  oldPosition[0][0] = initialPosition.getX();
  oldPosition[1][0] = initialPosition.getY();
  oldPosition[2][0] = initialPosition.getZ();
  
//   PrintMatrix(oldPosition);
  
  std::vector<std::vector<double> > newPosition;
  newPosition = MatrixMultiplication(transitionMatrix,oldPosition);
  
  G4ThreeVector positionAfterRotation;
  positionAfterRotation.setX(newPosition[0][0]);
  positionAfterRotation.setY(newPosition[1][0]);
  positionAfterRotation.setZ(newPosition[2][0]);
  
  return positionAfterRotation;

}


G4ThreeVector TruebeamHead::GetCenterPositionOfRotatedMLC(G4ThreeVector initialPosition, G4double gantryRtn, G4double collRtn)
{

  std::vector<std::vector<double> > transitionMatrixRg; // transition matrix of gantry angle
  std::vector<std::vector<double> > transitionMatrixRc; // transition matrix of collimation angle
  transitionMatrixRg =   MatrixInitialization(3,3);
  transitionMatrixRc =   MatrixInitialization(3,3);
    
//     transitionMatrixRg[0][0] = 1;
//     transitionMatrixRg[0][1] = 0;
//     transitionMatrixRg[0][2] = 0;
//     transitionMatrixRg[1][0] = 0;
//     transitionMatrixRg[1][1] = cos(gantryRtn);
//     transitionMatrixRg[1][2] = -sin(gantryRtn);
//     transitionMatrixRg[2][0] = 0;
//     transitionMatrixRg[2][1] = sin(gantryRtn);
//     transitionMatrixRg[2][2] = cos(gantryRtn); // gantry rotates with X
    
        
    transitionMatrixRg[0][0] = cos(gantryRtn);
    transitionMatrixRg[0][1] = 0;
    transitionMatrixRg[0][2] = sin(gantryRtn);
    transitionMatrixRg[1][0] = 0;
    transitionMatrixRg[1][1] = 1;
    transitionMatrixRg[1][2] = 0;
    transitionMatrixRg[2][0] = -sin(gantryRtn);
    transitionMatrixRg[2][1] = 0;
    transitionMatrixRg[2][2] = cos(gantryRtn); // gantry rotates with Y
    
    
    
    transitionMatrixRc[0][0] = cos(collRtn);
    transitionMatrixRc[0][1] = -sin(collRtn);
    transitionMatrixRc[0][2] = 0;
    transitionMatrixRc[1][0] = sin(collRtn);
    transitionMatrixRc[1][1] = cos(collRtn);
    transitionMatrixRc[1][2] = 0;
    transitionMatrixRc[2][0] = 0;
    transitionMatrixRc[2][1] = 0;
    transitionMatrixRc[2][2] = 1; // collimation rotates with Z
    
  
//   PrintMatrix(transitionMatrix);
  
    std::vector<std::vector<double> > transitionMatrixRg0; // transition matrix of gantry angle
    std::vector<std::vector<double> > transitionMatrixRc0; // transition matrix of collimation angle
    transitionMatrixRg0 =   MatrixInitialization(3,3);
    transitionMatrixRc0 =   MatrixInitialization(3,3);
    G4double gantryRtn0 = -pi;
//     transitionMatrixRg0[0][0] = 1;
//     transitionMatrixRg0[0][1] = 0;
//     transitionMatrixRg0[0][2] = 0;
//     transitionMatrixRg0[1][0] = 0;
//     transitionMatrixRg0[1][1] = cos(gantryRtn0);
//     transitionMatrixRg0[1][2] = -sin(gantryRtn0);
//     transitionMatrixRg0[2][0] = 0;
//     transitionMatrixRg0[2][1] = sin(gantryRtn0);
//     transitionMatrixRg0[2][2] = cos(gantryRtn0);
    
    
    
    transitionMatrixRg0[0][0] = cos(gantryRtn0);
    transitionMatrixRg0[0][1] = 0;
    transitionMatrixRg0[0][2] = sin(gantryRtn0);
    transitionMatrixRg0[1][0] = 0;
    transitionMatrixRg0[1][1] = 1;
    transitionMatrixRg0[1][2] = 0;
    transitionMatrixRg0[2][0] = -sin(gantryRtn0);
    transitionMatrixRg0[2][1] = 0;
    transitionMatrixRg0[2][2] = cos(gantryRtn0); // gantry rotates with Y
    
    G4double collRtn0 = -pi;
    
    transitionMatrixRc0[0][0] = cos(collRtn0);
    transitionMatrixRc0[0][1] = -sin(collRtn0);
    transitionMatrixRc0[0][2] = 0;
    transitionMatrixRc0[1][0] = sin(collRtn0);
    transitionMatrixRc0[1][1] = cos(collRtn0);
    transitionMatrixRc0[1][2] = 0;
    transitionMatrixRc0[2][0] = 0;
    transitionMatrixRc0[2][1] = 0;
    transitionMatrixRc0[2][2] = 1;
    
  
  std::vector<std::vector<double> > oldPosition;
  oldPosition = MatrixInitialization(3,1);
  oldPosition[0][0] = initialPosition.getX();
  oldPosition[1][0] = initialPosition.getY();
  oldPosition[2][0] = initialPosition.getZ();
  
//   PrintMatrix(oldPosition);
  
  
  
//   oldPosition = MatrixMultiplication(transitionMatrixRc0,oldPosition); // restore to original position
//   oldPosition = MatrixMultiplication(transitionMatrixRg0,oldPosition); // restore to original position
  
  std::vector<std::vector<double> > newPosition;
  newPosition = MatrixMultiplication(transitionMatrixRc,oldPosition); //implement collimation rotation
  newPosition = MatrixMultiplication(transitionMatrixRg,newPosition); // implement gantry rotation rotation
  
  
  G4ThreeVector positionAfterRotation;
  positionAfterRotation.setX(newPosition[0][0]);
  positionAfterRotation.setY(newPosition[1][0]);
  positionAfterRotation.setZ(newPosition[2][0]);
  
  return positionAfterRotation;
}

G4double TruebeamHead::MLCPosCorrection(G4double x)
{
  
  G4double p1 =  0.008 ;
  G4double p2 = 0;
  G4double p3 = 0;
  
//   return -1*(p1*x*x + p2*x + p3)/2;
  
  return 0;

}


void TruebeamHead::MLC()
{
    G4double MLCLength = 140*mm;
    G4double MLCEndCircleRadius = 234*mm;
    G4double MLCHeight = 77*mm;
    G4double theta = asin(MLCHeight/2/MLCEndCircleRadius);
    G4double MLCEndCircleAngle = 2*theta;
         
    G4double topSurfaceProximalMLCToISO  = 72.11*cm;
    G4double z_ProximalMLC = topSurfaceProximalMLCToISO - MLCHeight/2;
  
    G4double topSurfaceDistalMLCToISO = 63.39*cm;
    G4double z_DistalMLC = topSurfaceDistalMLCToISO - MLCHeight/2;
    
    z_ProximalMLC = 50*cm;
    
//     z_ProximalMLC = z_ProximalMLC - 5*cm; // this is for debugging the mlc position influence
    
    
    z_DistalMLC = z_ProximalMLC - 1*mm - MLCHeight;
    
    G4double MLCWidth = 1*cm;
    
    G4double MLCWidth_P = (SSD - (z_ProximalMLC- MLCHeight/4.0))/SSD*MLCWidth;
    
    G4double MLCWidth_D = (SSD - (z_DistalMLC - MLCHeight/4.0))/SSD*MLCWidth; // convert projected width to width at certain height
    
    cout<<"MLCWidth = "<<MLCWidth<<endl;
    cout<<"MLCWidth_P = "<<MLCWidth_P<<endl;
    cout<<"MLCWidth_D = "<<MLCWidth_D<<endl;
//     exit(0);
    
    G4Tubs* tempTubs_P = new G4Tubs("tempTubs", 0, MLCEndCircleRadius, 0.5*MLCWidth_P, pi/2- theta, MLCEndCircleAngle);
    G4Tubs* tempTubs_D = new G4Tubs("tempTubs", 0, MLCEndCircleRadius, 0.5*MLCWidth_D, pi/2- theta, MLCEndCircleAngle);
    
/*    
    G4GenericTrap(const G4String& pName,
    G4double pDz,
    const std::vector<G4TwoVector>& vertices)*/
  G4double pDz = MLCWidth;
  G4double x1 = 0;
  G4double y1 = 0;

  G4double x2 = MLCEndCircleRadius*cos(pi/2 + theta);
  G4double y2 =  MLCEndCircleRadius*sin(pi/2 + theta);
  G4double x3 = MLCEndCircleRadius*cos(pi/2 - theta);
  G4double y3 = MLCEndCircleRadius*sin(pi/2 - theta);
  
  cout<<y2<<endl;
//   exit(0);
  G4TwoVector p1 = G4TwoVector(x1,y1);
  G4TwoVector p2 = G4TwoVector(x2,y2);
  G4TwoVector p3 = G4TwoVector(x3,y3);
  G4TwoVector p4 = p1;
  G4TwoVector p5 = p1;
  G4TwoVector p6 = p2;
  G4TwoVector p7 = p3;
  G4TwoVector p8 = p4;
  
  
  std::vector<G4TwoVector> vertices;
  vertices.push_back(p1);
  vertices.push_back(p2);
  vertices.push_back(p3);
  vertices.push_back(p4);
  vertices.push_back(p5);
  vertices.push_back(p6);
  vertices.push_back(p7);
  vertices.push_back(p8);
  
  G4GenericTrap* tempTrap = new  G4GenericTrap("tempTrap",pDz,vertices);
    
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* tungstenAlloy = nist->FindOrBuildMaterial("tungstenAlloy");
  
  
  G4LogicalVolume*	 tempTubs_P_log = new G4LogicalVolume(tempTubs_P,tungstenAlloy,"tempTubs_P_log");
  G4LogicalVolume*	 tempTubs_D_log = new G4LogicalVolume(tempTubs_D,tungstenAlloy,"tempTubs_D_log");
  G4LogicalVolume* tempTrap_log = new G4LogicalVolume(tempTrap,tungstenAlloy,"tempTrap_log");
  
  
//   G4VSolid* substract = new G4SubtractionSolid("tempTubs-tempTrap",tempTubs,tempTrap,0,G4ThreeVector(0,0,0)); // this will have some geometry variations, not qutie right way
//   G4LogicalVolume* substract_log = new G4LogicalVolume(substract, tungstenAlloy,"substract_log");
  G4Box* tempBox_P = new G4Box("tempBox",MLCHeight/2, y2/2, MLCWidth_P/2);
  G4Box* tempBox_D = new G4Box("tempBox",MLCHeight/2, y2/2, MLCWidth_D/2);
  G4Transform3D transForUnio = G4Transform3D(G4RotationMatrix(),G4ThreeVector(0,y2/2, 0));
  
  G4VSolid* MLCLeaf_P = new G4UnionSolid("tempTubs + tempBox", tempTubs_P,tempBox_P, transForUnio);
  G4VSolid* MLCLeaf_D = new G4UnionSolid("tempTubs + tempBox", tempTubs_D,tempBox_D, transForUnio);
  //   G4VSolid* MLCLeaf = new G4Tubs("tempTubs", 0, MLCEndCircleRadius, 0.5*MLCWidth, 0, 2*pi);
  
  G4LogicalVolume* MLCLeaf_logPA = new G4LogicalVolume(MLCLeaf_P,tungstenAlloy,"MLCLeaf_log");  
    
  G4LogicalVolume* MLCLeaf_logPB = new G4LogicalVolume(MLCLeaf_P,tungstenAlloy,"MLCLeaf_log");
  G4LogicalVolume* MLCLeaf_logDA = new G4LogicalVolume(MLCLeaf_D,tungstenAlloy,"MLCLeaf_log");
  G4LogicalVolume* MLCLeaf_logDB = new G4LogicalVolume(MLCLeaf_D,tungstenAlloy,"MLCLeaf_log");
  
  
    
    G4RotationMatrix rotmA = G4RotationMatrix();
    rotmA.rotateZ(-90*deg);
    rotmA.rotateX(90*deg);
//     
    G4RotationMatrix rotmB = G4RotationMatrix();
    rotmB.rotateZ(90*deg);
    rotmB.rotateX(90*deg);
//     rotmB.rotateY(90*deg);
    G4double deltaX = MLCEndCircleRadius;
    
    
    
    ////////////////////////////////////////////////////////////////////
    /////////Read the Lianc information File to set up the MLC
    
     G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
     ReadLinacRTInformation readHlcy;
     readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information

    G4double gantryRtn = readHlcy.GetGantryRtnInformation()*deg;
    G4double collRtn = readHlcy.GetCollimationRtnInformation()*deg + 0*deg;
      
     cout<<"GantryRtn: "<<gantryRtn<<endl;
      cout<<"CollRtn: " << collRtn<<endl;
      
//       gantryRtn = gantryRtn - 180*deg;
//       collRtn = collRtn - 180*deg;
//       gantryRtn = -90*deg;
    std::vector<G4double> proximalDistVec;
    std::vector<G4double> distalDistVec;
    
    proximalDistVec = readHlcy.GetMLCProxInformation();
    distalDistVec = readHlcy.GetMLCDistInformation();
    
    for (int i =0; i<proximalDistVec.size(); i++)
    {
      proximalDistVec[i] = proximalDistVec[i]*cm; // set units 
    }
    
    for (int i = 0;i<distalDistVec.size(); i++)
    {
      distalDistVec[i] = distalDistVec[i]*cm;
    }
    
     cout<<"This is the mlc position"<<endl   ;
    for (int i = 0; i<proximalDistVec.size(); i++)
    {
      cout<<proximalDistVec[i]/cm<<endl;
    }
    
    
    for (int i = 0; i<distalDistVec.size(); i++)
    {
      cout<<distalDistVec[i]/cm <<endl;
    }
    
    
    // here we add some correction for the mlc position based on Varian's data
    
//     for (int i =0; i<proximalDistVec.size(); i++)
//     {
//       cout<<MLCPosCorrection(proximalDistVec[i]/cm)*mm<<endl;
//     }
//     
//     for (int i = 0;i<distalDistVec.size(); i++)
//     {
//       cout<<MLCPosCorrection(distalDistVec[i]/cm)*mm<<endl;
//     }
    

    for (int i =0; i<proximalDistVec.size(); i++)
    {
      proximalDistVec[i] = proximalDistVec[i] - MLCPosCorrection(proximalDistVec[i]/cm)*mm;
    }
    
    for (int i = 0;i<distalDistVec.size(); i++)
    {
      distalDistVec[i] = distalDistVec[i] - MLCPosCorrection(distalDistVec[i]/cm)*mm;
    }
    
    
    
    
    
    
//         // placing the proximal MLC , directly reading the MLC position file
//     
//     ifstream proxMLCFile("Prox.txt");
//     ifstream distMLCFile("Dist.txt");
//     
//     
//     stringstream proximalDistance;
//     stringstream distalDistance;
//     string line;
//     if (proxMLCFile.is_open())
//     {
//       
//       while (getline(proxMLCFile,line))
//       {
// // 	cout<<" proximal: "<<line<<endl;
// 	proximalDistance<<line;
//       }
//     }
//     string distance;
// 
//     while(getline(proximalDistance,distance,' '))
//     {
//       proximalDistVec.push_back(atof(distance.c_str())*cm);
//     }
//     
// 
//     if (distMLCFile.is_open())
//     {
//       while (getline(distMLCFile,line))
//       {
// // 	cout<<"distal: "<<line<<endl;
// 	distalDistance<<line;
//       }
//     }
//     
//     while(getline(distalDistance,distance,' '))
//     {
//       distalDistVec.push_back(atof(distance.c_str())*cm);
//     }
//     
    cout<<"This is the mlc position after correction"<<endl;
    for (int i = 0; i<proximalDistVec.size(); i++)
    {
      cout<<proximalDistVec[i]/cm<<endl;
    }
    
    
    for (int i = 0; i<distalDistVec.size(); i++)
    {
      cout<<distalDistVec[i]/cm<<endl;
    }

    

    cout<<"number of mlc leaves of proximal is "<<proximalDistVec.size()<<endl;
    cout<<"number of mlc leaves of distal is "<<distalDistVec.size()<<endl;
//     exit(0);


    
    
    for (int i = 0; i<29; i++)
    {
//       cout<<"proximal "<<proximalDistVec[i]<<endl;
      proximalDistVec[i] = GetMLCPositionA(proximalDistVec[i],SSD,MLCEndCircleRadius,z_ProximalMLC);
//       cout<<"proximial modifying "<<proximalDistVec[i]/mm<<endl;
    }
    
    for (int i = 0; i<29; i++)
    {
//       cout<<"proximal "<<proximalDistVec[i+29]<<endl;
      proximalDistVec[i+29] = GetMLCPositionB(proximalDistVec[i+29],SSD,MLCEndCircleRadius,z_ProximalMLC);
//       cout<<"proximial modifying "<<proximalDistVec[i+29]/mm<<endl;
    }
    
    
    for (int i = 0; i<28; i++)
    {
      distalDistVec[i] = GetMLCPositionA(distalDistVec[i],SSD,MLCEndCircleRadius,z_DistalMLC);
//       cout<<"distal "<<distalDistVec[i]/mm<<endl;
    }
    
        
    for (int i = 0; i<28; i++)
    {
      distalDistVec[i+28] = GetMLCPositionB(distalDistVec[i+28],SSD,MLCEndCircleRadius,z_DistalMLC);
//       cout<<"distal "<<distalDistVec[i+28]/mm<<endl;
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////Place the MLC leaves in to system


//     for (int i = 0; i<1; i++)
    for (int i = 0; i<29; i++)
    {
      G4ThreeVector oldPosition = G4ThreeVector(proximalDistVec[28-i] ,0 + 14*MLCWidth_P - i*MLCWidth_P ,z_ProximalMLC);
      G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn, collRtn);
      G4RotationMatrix rotmAtmp = rotmA;
//       rotmAtmp.rotateY(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
//       rotmAtmp.rotateZ(-180*deg);// Get back to the orginal rotation from current angle
//       rotmAtmp.rotateY(gantryRtn);
//       rotmAtmp.rotateZ(collRtn);
	rotmAtmp.rotateZ(collRtn);
	rotmAtmp.rotateY(gantryRtn);
      G4Transform3D transform = G4Transform3D(rotmAtmp,newPosition);
//       G4ThreeVector invertNewPosition = G4ThreeVector(-newPosition.getX(),newPosition.getY(),newPosition.getZ());
//       G4Transform3D transform = G4Transform3D(rotmAtmp,invertNewPosition);
      new G4PVPlacement(transform,"unionV", MLCLeaf_logPA,PVWorld,false,0,false); // place the proximal A
    }
//     
//        
//     for (int i = 0; i<1; i++)
    for (int i = 0; i<29; i++)
    {
      G4ThreeVector oldPosition = G4ThreeVector(proximalDistVec[57-i],0 + 14*MLCWidth_P - i*MLCWidth_P ,z_ProximalMLC);
      G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn,collRtn);
      G4RotationMatrix rotmBtmp = rotmB;
//       rotmBtmp.rotateY(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
//       rotmBtmp.rotateZ(-180*deg);
//       rotmBtmp.rotateY(-180*deg);
//       rotmBtmp.rotateY(gantryRtn);
//       rotmBtmp.rotateZ(collRtn);
	rotmBtmp.rotateZ(collRtn);
	rotmBtmp.rotateY(gantryRtn);
	
     G4Transform3D transform = G4Transform3D(rotmBtmp,newPosition);
//       G4ThreeVector invertNewPosition = G4ThreeVector(-newPosition.getX(),newPosition.getY(),newPosition.getZ());
//       G4Transform3D transform = G4Transform3D(rotmBtmp,invertNewPosition);
      new G4PVPlacement(transform,"unionV", MLCLeaf_logPB,PVWorld,false,0,false); // place the proximal B
    }

       
   for (int i = 0; i<28; i++)
   {
      G4ThreeVector oldPosition = G4ThreeVector(distalDistVec[27-i],0 + MLCWidth_D/2 + 13*MLCWidth_D - i*MLCWidth_D ,z_DistalMLC );
      G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn,collRtn); 
      G4RotationMatrix rotmAtmp = rotmA;
//       rotmAtmp.rotateY(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
//       rotmAtmp.rotateZ(-180*deg);// Get back to the orginal rotation from current angle
//       rotmAtmp.rotateY(gantryRtn);
//       rotmAtmp.rotateZ(collRtn);
	rotmAtmp.rotateZ(collRtn);
	rotmAtmp.rotateY(gantryRtn);
     
      
      G4Transform3D transform = G4Transform3D(rotmAtmp,newPosition);
      new G4PVPlacement(transform,"unionV", MLCLeaf_logDA,PVWorld,false,0,false); // place the distal A
   }

   for (int i = 0; i<28; i++)
   {
     G4ThreeVector oldPosition = G4ThreeVector(distalDistVec[55-i],0 + MLCWidth_D/2 + 13*MLCWidth_D - i*MLCWidth_D ,z_DistalMLC );
     G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn,collRtn); 
      G4RotationMatrix rotmBtmp = rotmB;
//       rotmBtmp.rotateY(-180*deg);// Get back to the orginal rotation from current angle, current is 180 deg, 180 deg
//       rotmBtmp.rotateZ(-180*deg);
//       rotmBtmp.rotateY(-180*deg);
//       rotmBtmp.rotateY(gantryRtn);
//       rotmBtmp.rotateZ(collRtn);
        rotmBtmp.rotateZ(collRtn);
	rotmBtmp.rotateY(gantryRtn);
    
      
     G4Transform3D transform = G4Transform3D(rotmBtmp,newPosition);
      new G4PVPlacement(transform,"unionV", MLCLeaf_logDB,PVWorld,false,0,false);  // place the distal B
   }


///////////////////////////////////////////////////////////////////////////////
////////////// old design, fixed gantry rotation, and collimation rotatioin, pi, pi

    
//     for (int i = 0; i<29; i++)
//     {
// //       G4Transform3D transform = G4Transform3D(rotmA,G4ThreeVector(0 - deltaX  - proximalDistVec[i] ,0 + 14*MLCWidth - i*MLCWidth ,z_ProximalMLC));
//       G4Transform3D transform = G4Transform3D(rotmA,G4ThreeVector(proximalDistVec[i] ,0 + 14*MLCWidth - i*MLCWidth ,z_ProximalMLC));
//       new G4PVPlacement(transform,"unionV", MLCLeaf_logPA,PVWorld,false,0,false); // place the proximal A
//     }
//     
//        
//     for (int i = 0; i<29; i++)
//     {
//       G4Transform3D transform = G4Transform3D(rotmB,G4ThreeVector(proximalDistVec[i+29],0 + 14*MLCWidth - i*MLCWidth ,z_ProximalMLC));
//       new G4PVPlacement(transform,"unionV", MLCLeaf_logPB,PVWorld,false,0,false); // place the proximal B
//     }
// 
//        
//    for (int i = 0; i<28; i++)
//    {
//       G4Transform3D transform = G4Transform3D(rotmA,G4ThreeVector(distalDistVec[i],0 + MLCWidth/2 + 13*MLCWidth - i*MLCWidth ,z_DistalMLC ));
//       new G4PVPlacement(transform,"unionV", MLCLeaf_logDA,PVWorld,false,0,false); // place the distal A
//    }
// 
//    for (int i = 0; i<28; i++)
//    {
//       G4Transform3D transform = G4Transform3D(rotmB,G4ThreeVector(distalDistVec[i+28],0 + MLCWidth/2 + 13*MLCWidth - i*MLCWidth ,z_DistalMLC ));
//       new G4PVPlacement(transform,"unionV", MLCLeaf_logDB,PVWorld,false,0,false);  // place the distal B
//    }





        // Region for cuts
    G4Region *regVol;
    regVol= new G4Region("MLCRegion");
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(2.*cm); // set the production cut of MLC
    regVol->SetProductionCuts(cuts);
    MLCLeaf_logPA->SetRegion(regVol);
    MLCLeaf_logPB->SetRegion(regVol);
    MLCLeaf_logDA->SetRegion(regVol);
    MLCLeaf_logDB->SetRegion(regVol);
    
//     MLCLeaf_logPA->SetUserLimits(new G4UserLimits(20*cm)); // here we set up the maximum step size
//     MLCLeaf_logPB->SetUserLimits(new G4UserLimits(20*cm));
//     MLCLeaf_logDA->SetUserLimits(new G4UserLimits(20*cm));
//     MLCLeaf_logDB->SetUserLimits(new G4UserLimits(20*cm));
    
    G4VisAttributes*  simpleAlSVisAttPA= new G4VisAttributes(G4Colour::Blue());
    simpleAlSVisAttPA->SetVisibility(true);
    simpleAlSVisAttPA->SetForceSolid(true);
    MLCLeaf_logPA->SetVisAttributes(simpleAlSVisAttPA);
    
    G4VisAttributes*  simpleAlSVisAttPB= new G4VisAttributes(G4Colour::Cyan());
    simpleAlSVisAttPB->SetVisibility(true);
    simpleAlSVisAttPB->SetForceSolid(true);
    MLCLeaf_logPB->SetVisAttributes(simpleAlSVisAttPB);
    
    G4VisAttributes*  simpleAlSVisAttDA= new G4VisAttributes(G4Colour::Red());
    simpleAlSVisAttDA->SetVisibility(true);
    simpleAlSVisAttDA->SetForceSolid(true);
    MLCLeaf_logDA->SetVisAttributes(simpleAlSVisAttDA);
    
    G4VisAttributes*  simpleAlSVisAttDB= new G4VisAttributes(G4Colour::Magenta());
    simpleAlSVisAttDB->SetVisibility(true);
    simpleAlSVisAttDB->SetForceSolid(true);
    MLCLeaf_logDB->SetVisAttributes(simpleAlSVisAttDB);

}

G4double TruebeamHead::GetMLCPositionA(G4double kx,G4double ssd, G4double R,  G4double y)
{
  G4double L = -kx; // kx is the projected MLC position
  G4double k = -ssd/L;
  cout<<"A bank K = "<<k<<endl;
  G4double x ;
  
  x = (y - sqrt(1+ k*k)*R)/k + L;
  if (kx>0)
  {
     return x;
  }
  else
  {
//     x = (y + sqrt(1+ k*k)*R)/k + L;
    return x - 2*R;
  }
  
//   return x;

}
G4double TruebeamHead::GetMLCPositionB(G4double kx, G4double ssd, G4double R, G4double y)
{
  G4double L = kx;
  G4double k = -ssd/L;
  cout<<"B bank k = "<<k<<endl;
  G4double x;
   x = (y-sqrt(1+ k*k)*R)/k + L;
  if (kx>0)
  {
   return x;
  }
  else
  {
//     x = (y+sqrt(1+ k*k)*R)/k + L;
    return x + 2*R;
  }
//   return x;

}




std::vector< std::vector< double > > TruebeamHead::MatrixInitialization(int m, int n)
{
  
    std::vector<std::vector<double> > M;
    for (int i=0;i<m;i++ )
    {
        M.push_back(std::vector<double>());
        for (int j=0;j<n;j++)
        {
            M[i].push_back(0); // initial result
        }
    }
    return M;
}
std::vector< std::vector< double > > TruebeamHead::MatrixMultiplication(std::vector< std::vector< double > > M1, std::vector< std::vector< double > > M2)
{
    std::vector<std::vector<double> >  result;
    int m = M1.size();
    int n = M2[0].size();
    result = MatrixInitialization(m,n);
//     cout<<"m = "<<m<<" n = "<<n<<endl;
//     PrintMatrix(result);
    
    for (int i = 0;i<M1.size();i++)
    {
        for (int j = 0;j< M2[0].size();j++)
        {
            for (int k=0;k<M1[0].size();k++)
            {
                result[i][j] += M1[i][k] * M2[k][j];
            }
        }
    }
    return result;  
}

void TruebeamHead::PrintMatrix(std::vector< std::vector< double > > M)
{
    for (int i =0;i<M.size();i++)
    {
        for (int j=0;j<M[i].size();j++)
        {
                cout<<"M"<<i<<","<<j<<" = "<<M[i][j]<<endl;
        }
    }

}


void TruebeamHead::SetFieldSize(G4double fs)
{
  fieldSize  = fs;

}
void TruebeamHead::SetSSD(G4double ssd)
{
  SSD = ssd ;
  
  targetPosition.set(0,0,ssd);

}

