#include "MRLinacHead.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <stdlib.h>
using namespace std;


MRLinacHead::MRLinacHead()
{


}

MRLinacHead::~MRLinacHead()
{

}



void MRLinacHead::Target()
{

  //    materials  

	G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material* Re = G4NistManager::Instance()->FindOrBuildMaterial("G4_Re");
	
	G4Material* WReAlloy = new G4Material("WReAlloy",20.54*g/cm3, 2); // define the Tungsten and Rhenium alloy material
	WReAlloy->AddMaterial(W,75.0*perCent);
	WReAlloy->AddMaterial(Re,25.0*perCent);
	

	//    colors

   G4Colour  cyan    (0.0, 1.0, 1.0);
   G4Colour  magenta (1.0, 0.0, 1.0); 
 
  //    volumes
  //    beam line along z axis
//------------------------target 6MV------------------------

  
  G4double targetTungstenThickness = 0.899*mm;
  G4double targetTungstenRadius = 7.0/2*mm;
  G4double targetCuThickness = 2.8*mm;
  G4Tubs* targetTungstenTube = new G4Tubs("targetTungsten",0.,targetTungstenRadius,0.5*targetTungstenThickness,0.,twopi);
  
  G4LogicalVolume *targetTungsten_log = new G4LogicalVolume(targetTungstenTube,WReAlloy,"targetTungsten_log",0,0,0);
  
  G4double targetTungstenPos_x = 0.0*mm;
  G4double targetTungstenPos_y = 0.0*mm;
  G4double beamPortIsoDistance  = 90*cm;
  G4double targetTungstenPos_z = (beamPortIsoDistance + targetCuThickness+0.5*targetTungstenThickness);
  targetTungsten_phys = new G4PVPlacement(0,
            G4ThreeVector(targetTungstenPos_x,targetTungstenPos_y,targetTungstenPos_z),
            "targetTungsten",targetTungsten_log,PVWorld,false,0);

  
  
  G4Tubs* targetCuTube = new G4Tubs("targetCu",0.,targetTungstenRadius,0.5*targetCuThickness,0.,twopi);
  
  G4LogicalVolume *targetCu_log = new G4LogicalVolume(targetCuTube,Cu,"targetCu_log",0,0,0);
  G4double targetCuPos_x = 0.0*m;
  G4double targetCuPos_y = 0.0*m;
  G4double targetCuPos_z = targetTungstenPos_z - 0.5*(targetTungstenThickness+targetCuThickness);
  
  targetPosition.setX(targetCuPos_x);
  targetPosition.setY(targetCuPos_y);
  targetPosition.setZ(targetCuPos_z);
  
  targetCu_phys = new G4PVPlacement(0,
            G4ThreeVector(targetCuPos_x,targetCuPos_y,targetCuPos_z),
            "targetCu",targetCu_log,PVWorld,false,0);


	// ***********  REGIONS for CUTS

  	G4Region *regVol;
	regVol= new G4Region("targetR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	targetTungsten_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(targetTungsten_log);
	targetCu_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(targetCu_log);

	// Visualization attributes

	G4VisAttributes* simpleWSVisAtt, *simpleCuSVisAtt;
	simpleWSVisAtt= new G4VisAttributes(magenta);
	simpleWSVisAtt->SetVisibility(true);
// 	simpleWSVisAtt->SetForceSolid(true);
	simpleCuSVisAtt= new G4VisAttributes(cyan);
	simpleCuSVisAtt->SetVisibility(true);
// 	simpleCuSVisAtt->SetForceSolid(true);
	targetTungsten_log->SetVisAttributes(simpleWSVisAtt);
	targetCu_log->SetVisAttributes(simpleCuSVisAtt);

}

void MRLinacHead::PrimaryCollimator()
{
  //    materials 

	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4Material* Ni = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
	G4Material* default_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"); // take it as vacuum
	G4Material* WCuNiAlloy = new G4Material("WCuNiAlloy",17.75*g/cm3,3);// define the Tungsten Cu Nickel alloy
	G4double fw = 75.0*perCent; // suppose the Tungsten mass fraction is 75%
	G4double fCu = 15.57/2*perCent; // this is calcualted assuming mass fraction of Cu and Ni are equal
	G4double fNi = 1-fCu-fw;
	
	WCuNiAlloy->AddMaterial(W,fw);
	WCuNiAlloy->AddMaterial(Cu,fCu);
	WCuNiAlloy->AddMaterial(Ni,fNi);
	
	G4double primaryCollimatorBeamConeThickness = 101*mm;
	G4double baseRadius = (diameterOfPrimaryCollimatorTop  + primaryCollimatorBeamConeThickness*tan(maxConeAngle/2)*2)/2;
	G4double primaryCollimatorRadius = 4*baseRadius;
	G4double primaryCollimatorThickness = 101*mm;
	G4Tubs* primaryCollimatorTube = new G4Tubs("primaryCollimatorTube",0.,primaryCollimatorRadius,0.5*primaryCollimatorThickness,0.,twopi);
	G4LogicalVolume*  primaryCollimatorLogicV = new G4LogicalVolume(primaryCollimatorTube,WCuNiAlloy,"primaryCollimatorLogicV",0,0,0);
	
	G4double primaryCollimatorPos_x = 0.0*cm;
	G4double primaryCollimatorPos_y = 0.0*cm;
	G4double primaryCollimatorPordIsoDistance = 799*mm;
	G4double primaryCollimatorPos_z = (primaryCollimatorPordIsoDistance + 0.5*primaryCollimatorThickness);
	
	G4double primaryCollimatorBeamConeRmin1 = 0.*mm;
	G4double primaryCollimatorBeamConeRmax1 =  diameterOfPrimaryCollimatorTop/2;

	G4double primaryCollimatorBeamConeRmin2 = 0.0*mm;
	G4double primaryCollimatorBeamConeRmax2 = (diameterOfPrimaryCollimatorTop  + primaryCollimatorBeamConeThickness*tan(maxConeAngle/2)*2)/2;
	
// 	cout<<"The primaryCollimatorBeamConeRmin1 is "<<primaryCollimatorBeamConeRmax2<<endl;
// 	exit(0);
	
	
	G4Cons* primaryCollimatorBeamCone = new  G4Cons("primaryCollimatorBeamCone", primaryCollimatorBeamConeRmin2,primaryCollimatorBeamConeRmax2,
	  primaryCollimatorBeamConeRmin1,primaryCollimatorBeamConeRmax1,
	  primaryCollimatorBeamConeThickness/2,0.,twopi
	);
	
	G4VSolid* subtract = new G4SubtractionSolid("primaryCollimatorTube-primaryCollimatorBeamCone",primaryCollimatorTube,
	  primaryCollimatorBeamCone,0,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm)
	);
	
	G4LogicalVolume* primaryCollimatorWithConeLogicV = new G4LogicalVolume(subtract,WCuNiAlloy,"primaryCollimatorWithCone",0,0,0);
	
// 	primaryCollimator_phys = new G4PVPlacement(0,G4ThreeVector(primaryCollimatorPos_x,primaryCollimatorPos_y,primaryCollimatorPos_z),
// 	  "primaryCollimator",primaryCollimatorLogicV,PVWorld,false,0);
    
	
	primaryCollimator_phys = new G4PVPlacement(0,G4ThreeVector(primaryCollimatorPos_x,primaryCollimatorPos_y,primaryCollimatorPos_z),
	  "primaryCollimatorWithCone",primaryCollimatorWithConeLogicV,PVWorld,false,0);
    

	

  //    colors

   G4Colour  magenta (1.0, 0.0, 1.0); 
   

  
//--------- Visualization attributes -------------------------------
  G4VisAttributes* primaryCollimatorVisAtt = new G4VisAttributes(magenta);
  primaryCollimatorVisAtt->SetVisibility(true);
  primaryCollimatorWithConeLogicV->SetVisAttributes(primaryCollimatorVisAtt);
  
	// ***********  REGIONS for CUTS for PrimaryCollimator

  G4Region *regVol;
  regVol= new G4Region("PrymCollR");
  G4ProductionCuts* cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.1*cm);
  regVol->SetProductionCuts(cuts);

  primaryCollimatorWithConeLogicV->SetRegion(regVol);
  regVol->AddRootLogicalVolume(primaryCollimatorWithConeLogicV);
 

   

}

void MRLinacHead::ElectronAbsorber()
{
  G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4double electronAbsorberPordIsoDistance = 795*mm;
  G4double beamDistance = 900*mm-electronAbsorberPordIsoDistance;
  G4double electronAbsorberTubeThickness = 4*mm;
  G4double electronAbsorberTubeRadius = (diameterOfPrimaryCollimatorTop  + beamDistance*tan(maxConeAngle/2)*2)/2;
  G4Tubs* electronAbsorberTube = new G4Tubs("electronAbsorberTube",0.,electronAbsorberTubeRadius,0.5*electronAbsorberTubeThickness,0.,twopi);
  
  G4LogicalVolume* electronAbsorberLogicV = new G4LogicalVolume(electronAbsorberTube,Al,"electronAbsorberLogicV",0,0,0);
  
  G4double electronAbsorberPos_x = 0.0*cm;
  G4double electronAbsorberPos_y = 0.0*cm;
  G4double electronAbsorberPos_z = (electronAbsorberPordIsoDistance + 0.5*electronAbsorberTubeThickness);
  
  electronAbsorber_phys = new G4PVPlacement(0,G4ThreeVector(electronAbsorberPos_x,electronAbsorberPos_y,electronAbsorberPos_z),
    "electronAbsorber_phys",electronAbsorberLogicV,PVWorld,false,0
  );
  
    G4VisAttributes* simpleAlSVisAtt;
    simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
    simpleAlSVisAtt->SetVisibility(true);
    electronAbsorberLogicV->SetVisAttributes(simpleAlSVisAtt);



}

void MRLinacHead::IonizationChamber()
{
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4double ionizationChamberPordIsoDistance = 783*mm;
    G4double effectiveDistance = 900*mm- ionizationChamberPordIsoDistance;
    G4double ionizationChamberTubeRadius = (diameterOfPrimaryCollimatorTop  + effectiveDistance*tan(maxConeAngle/2)*2)/2;
    G4double ionizationChamberTubeThickness = 7.0*mm;
    G4Tubs* ionizationChamberTube = new G4Tubs("ionizationChamberTube",0.,ionizationChamberTubeRadius,
      0.5*ionizationChamberTubeThickness,0.,twopi
    );
    
    G4LogicalVolume* ionizationChamberLogicV = new G4LogicalVolume(ionizationChamberTube,Al,"ionizationChamberLogicV",0,0,0);
    G4double ionizationChamberPos_x = 0.0*mm;
    G4double ionizationChamberPos_y = 0.0*mm;

    G4double ionizationChamberPos_z = (ionizationChamberPordIsoDistance + 0.5*ionizationChamberTubeThickness);
    
    ionizationChamber_phys = new G4PVPlacement(0,G4ThreeVector(ionizationChamberPos_x,ionizationChamberPos_y,ionizationChamberPos_z),
      "ionizationChamber_phys",ionizationChamberLogicV,PVWorld,false,0
    );
    
    
    G4VisAttributes* simpleAlSVisAtt;
    simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
    simpleAlSVisAtt->SetVisibility(true);
    ionizationChamberLogicV->SetVisAttributes(simpleAlSVisAtt);

}


void MRLinacHead::MLC()
{
    //    materials 

	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4Material* Ni = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
	G4Material* default_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"); // take it as vacuum
	G4Material* WCuNiAlloy = new G4Material("WCuNiAlloy",17.75*g/cm3,3);// define the Tungsten Cu Nickel alloy
	G4double fw = 75.0*perCent; // suppose the Tungsten mass fraction is 75%
	G4double fCu = 15.57/2*perCent; // this is calcualted assuming mass fraction of Cu and Ni are equal
	G4double fNi = 1-fCu-fw;
	
	WCuNiAlloy->AddMaterial(W,fw);
	WCuNiAlloy->AddMaterial(Cu,fCu);
	WCuNiAlloy->AddMaterial(Ni,fNi);
	
	G4double firstMLCTopLayerPordIsoDistance = 635*mm;
	G4double beamPortIsoDistance = 900*mm;
	G4double firstMLCThickness = 55*mm;
	
	G4double firstMLCConeRmin2 = 0.*mm;
	G4double firstMLCConeRmax2 = (diameterOfPrimaryCollimatorTop  + (beamPortIsoDistance-firstMLCTopLayerPordIsoDistance)*tan(maxConeAngle/2)*2)/2;
	
// 	cout<<"firstMLCConeRmax2="<<firstMLCConeRmax2<<endl;
	
	
	G4double firstMLCConeRmin1 = 0.*mm;
	G4double firstMLCConeRmax1 = (diameterOfPrimaryCollimatorTop  + (beamPortIsoDistance-firstMLCTopLayerPordIsoDistance+ firstMLCThickness)*tan(maxConeAngle/2)*2)/2;
	
// 	cout<<"firstMLCConeRmax1="<<firstMLCConeRmax1<<endl;
// 	exit(0);
	G4Cons* firstMLCCone = new G4Cons("firstMLCCone",firstMLCConeRmin1,firstMLCConeRmax1,
	  firstMLCConeRmin2,firstMLCConeRmax2,firstMLCThickness/2,0.,twopi
	);
	
	G4double firstMLCConeWideRmin2 = 0.*mm;
	G4double firstMLCConeWideRmax2 = 2*firstMLCConeRmax2; // here we set the wide part is two times wider than the air coneAngle
	G4double firstMLCConeWideRmin1 = 0.*mm;
	G4double firstMLCConeWideRmax1 = 2*firstMLCConeRmax1;
	
	G4Cons* firstMLCConeWide= new G4Cons("firstMLCConeWide",firstMLCConeWideRmin1, firstMLCConeWideRmax1,
	  firstMLCConeWideRmin2, firstMLCConeWideRmax2,firstMLCThickness/2,0.,twopi
	);
	
	G4VSolid* subtract = new G4SubtractionSolid("firstMLCConeWide-firstMLCCone",firstMLCConeWide,firstMLCCone,
	  0,G4ThreeVector(0.*mm,0.*mm,0.*mm)
	);
	
	G4LogicalVolume* firstMLCConeLogicV = new G4LogicalVolume(subtract,WCuNiAlloy,"firstMLCConeLogicV",0,0,0);
	
	G4double firstMLCConePos_x = 0.*mm;
	G4double firstMLCConePos_y = 0.*mm;
	G4double firstMLCConePos_z = firstMLCTopLayerPordIsoDistance - 0.5*firstMLCThickness;
	firstMLC_phys = new G4PVPlacement(0,G4ThreeVector(firstMLCConePos_x,firstMLCConePos_y,firstMLCConePos_z),
	  "firstMLC_phys",firstMLCConeLogicV,PVWorld,false,0
	);
	

	G4double secondMLCTopLayerPordIsoDistance = firstMLCTopLayerPordIsoDistance - firstMLCThickness - 20*mm;
	G4double secondMLCConeRmin2 = 0.*mm;
	G4double secondMLCConeRmax2 = (diameterOfPrimaryCollimatorTop  + (beamPortIsoDistance-secondMLCTopLayerPordIsoDistance)*tan(maxConeAngle/2)*2)/2;
	
	G4double secondMLCConeRmin1 = 0.*mm;
	G4double secondMLCConeRmax1 = (diameterOfPrimaryCollimatorTop  + (beamPortIsoDistance-secondMLCTopLayerPordIsoDistance+ firstMLCThickness)*tan(maxConeAngle/2)*2)/2;
	
	G4Cons* secondMLCCOne = new G4Cons("secondMLCCone",secondMLCConeRmin1,secondMLCConeRmax1,
	  secondMLCConeRmin2,secondMLCConeRmax2,firstMLCThickness/2,0.,twopi
	);
	
	G4double secondMLCConeWideRmin2 = 0.*mm;
	G4double secondMLCConeWideRmax2 = 2*secondMLCConeRmax2;
	G4double secondMLCConeWideRMin1 = 0.*mm;
	G4double secondMLCConeWideRMax1 = 2*secondMLCConeRmax1;
	
	G4Cons* secondMLCConeWide = new G4Cons("secondMLCConeWide",secondMLCConeWideRMin1,secondMLCConeWideRMax1,
	  secondMLCConeWideRmin2,secondMLCConeWideRmax2,firstMLCThickness/2,0.,twopi
	);
	
	G4VSolid* subtract1 = new G4SubtractionSolid("secondMLCConeWide-secondMLCCone",secondMLCConeWide,secondMLCCOne,0,
	  G4ThreeVector(0.*mm,0.*mm,0.*mm)
	);
	
	G4LogicalVolume* secondMLCConeLogicV = new G4LogicalVolume(subtract1,WCuNiAlloy,"secondMLCConeLogicV",0,0,0);
	G4double secondMLCConePos_x = 0.*mm;
	G4double secondMLCConePos_y = 0.*mm;
	G4double secondMLCConePos_z = secondMLCTopLayerPordIsoDistance - 0.5*firstMLCThickness;
	secondMLC_phys = new G4PVPlacement(0,G4ThreeVector(secondMLCConePos_x,secondMLCConePos_y,secondMLCConePos_z),
	  "secondMLC_phys",secondMLCConeLogicV,PVWorld,false,0
	);
	
	
	

	
	
	
	
	G4VisAttributes* simpleAlSVisAtt;
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
	firstMLCConeLogicV->SetVisAttributes(simpleAlSVisAtt);
	secondMLCConeLogicV->SetVisAttributes(simpleAlSVisAtt);

}
void MRLinacHead::MLCSteepEdge()
{
  //    materials 

      G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
      G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
      G4Material* Ni = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
      G4Material* default_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"); // take it as vacuum
      G4Material* WCuNiAlloy = new G4Material("WCuNiAlloy",17.75*g/cm3,3);// define the Tungsten Cu Nickel alloy
      G4double fw = 75.0*perCent; // suppose the Tungsten mass fraction is 75%
      G4double fCu = 15.57/2*perCent; // this is calcualted assuming mass fraction of Cu and Ni are equal
      G4double fNi = 1-fCu-fw;
      
      WCuNiAlloy->AddMaterial(W,fw);
      WCuNiAlloy->AddMaterial(Cu,fCu);
      WCuNiAlloy->AddMaterial(Ni,fNi);
      

      G4double beamPortIsoDistance = 900*mm;
      G4double firstMLCThickness = 55*mm;
      

      G4double firstMLCInnerRadius = (diameterOfPrimaryCollimatorTop  + (beamPortIsoDistance)*tan(maxConeAngle/2)*2)/2;
      G4double firstMLCSolidTubeRadius = 2*firstMLCInnerRadius;
      
	cout<<"field size="<<firstMLCInnerRadius*2<<endl;
      
// 	exit(0);
      G4Tubs* firstMLC = new G4Tubs("firstMLC",firstMLCInnerRadius,firstMLCSolidTubeRadius,0.5*firstMLCThickness,0.,twopi);
      
      G4LogicalVolume* firstMLCLogicV = new G4LogicalVolume(firstMLC,WCuNiAlloy,"firstMLCLogicV",0,0,0);
      
      G4double firstMLCTopLayerPordIsoDistance = 635*mm;
      G4double firstMLCPos_x = 0.*mm;
      G4double firstMLCPos_y = 0.*mm;
      G4double firstMLCPos_z = firstMLCTopLayerPordIsoDistance - 0.5*firstMLCThickness;
      firstMLC_phys = new G4PVPlacement(0,G4ThreeVector(firstMLCPos_x,firstMLCPos_y,firstMLCPos_z),
	"firstMLC_phys",firstMLCLogicV,PVWorld,false,0
      );
      

      G4double secondMLCTopLayerPordIsoDistance = firstMLCTopLayerPordIsoDistance - firstMLCThickness - 20*mm;


      G4double secondMLCPos_x = 0.*mm;
      G4double secondMLCPos_y = 0.*mm;
      G4double secondMLCPos_z = secondMLCTopLayerPordIsoDistance - 0.5*firstMLCThickness;
      secondMLC_phys = new G4PVPlacement(0,G4ThreeVector(secondMLCPos_x,secondMLCPos_y,secondMLCPos_z),
	"secondMLC_phys",firstMLCLogicV,PVWorld,false,0
      );
      
      G4VisAttributes* simpleAlSVisAtt;
      simpleAlSVisAtt= new G4VisAttributes(G4Colour::Green());
      simpleAlSVisAtt->SetVisibility(true);
      firstMLCLogicV->SetVisAttributes(simpleAlSVisAtt);


}


void MRLinacHead::MLCCustomized(G4double coneAngle)
{
    //    materials 

      G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
      G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
      G4Material* Ni = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
      G4Material* default_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"); // take it as vacuum
      G4Material* WCuNiAlloy = new G4Material("WCuNiAlloy",17.75*g/cm3,3);// define the Tungsten Cu Nickel alloy
      G4double fw = 75.0*perCent; // suppose the Tungsten mass fraction is 75%
      G4double fCu = 15.57/2*perCent; // this is calcualted assuming mass fraction of Cu and Ni are equal
      G4double fNi = 1-fCu-fw;
      
      WCuNiAlloy->AddMaterial(W,fw);
      WCuNiAlloy->AddMaterial(Cu,fCu);
      WCuNiAlloy->AddMaterial(Ni,fNi);
      
      G4double referenceBeamPortIsoDistance = 900*mm + diameterOfPrimaryCollimatorTop/2/tan(maxConeAngle/2);
//       G4double beamPortIsoDistance = 900*mm;
      G4double beamPortIsoDistance = referenceBeamPortIsoDistance;
      G4double firstMLCTopLayerPordIsoDistance = 635*mm;
      G4double firstMLCThickness = 55*mm;
      

      G4double firstMLCInnerRadius =  beamPortIsoDistance - firstMLCTopLayerPordIsoDistance;
      G4double firstMLCOuterRadius = firstMLCInnerRadius+ firstMLCThickness;

      
// 	exit(0);
      G4double maximumFieldSizeDimX = 27*cm;
      G4double MLCArcAngle = twopi/10;
      G4Tubs* firstMLC = new G4Tubs("firstMLC",firstMLCInnerRadius,firstMLCOuterRadius,0.5*maximumFieldSizeDimX,0.,MLCArcAngle);
      
      G4LogicalVolume* firstMLCLogicV = new G4LogicalVolume(firstMLC,WCuNiAlloy,"firstMLCLogicV",0,0,0);
      
      

      G4double firstMLCPos_x = 0.*mm;
      G4double firstMLCPos_y = 0.*mm;
      G4double firstMLCPos_z = beamPortIsoDistance;
      G4RotationMatrix rotm  = G4RotationMatrix();
      rotm.rotateY(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis

      rotm.rotateX(coneAngle/2); // for first one, we let it rotate in clockwise
//       rotm.rotateX(0); // for first one, we let it rotate in clockwise
      
      G4ThreeVector position = G4ThreeVector(firstMLCPos_x,firstMLCPos_y,firstMLCPos_z);
      G4Transform3D transform = G4Transform3D(rotm,position);
      firstMLC_phys = new G4PVPlacement(transform,             //rotation,position
				  "firstMLC_phys",
				  firstMLCLogicV,
				  PVWorld,
				  false,
				  0,
				  false);
      
      rotm.rotateX(-MLCArcAngle - coneAngle); // for third one, we let it rotate anti-clockwise
      transform = G4Transform3D(rotm,position);
      thirdMLC_phys = new G4PVPlacement(transform,             //rotation,position
				  "firstMLC_phys",
				  firstMLCLogicV,
				  PVWorld,
				  false,
				  0,
				  false);


      G4double secondMLCTopLayerPordIsoDistance = firstMLCTopLayerPordIsoDistance - firstMLCThickness - 20*mm;
      
      G4double secondMLCInnerRadius = beamPortIsoDistance - secondMLCTopLayerPordIsoDistance;
      G4double secondMLCOuterRadius = secondMLCInnerRadius + firstMLCThickness;


      G4Tubs* secondMLC = new G4Tubs("firstMLC",secondMLCInnerRadius,secondMLCOuterRadius,0.5*maximumFieldSizeDimX,0.,MLCArcAngle);
      G4LogicalVolume* secondMLCLogicV = new G4LogicalVolume(secondMLC,WCuNiAlloy,"firstMLCLogicV",0,0,0);

      G4double secondMLCPos_x = 0.*mm;
      G4double secondMLCPos_y = 0.*mm;
      G4double secondMLCPos_z = beamPortIsoDistance;
      G4RotationMatrix rotm1  = G4RotationMatrix();
      rotm1.rotateY(90*deg);  // rotating from Y Z plane to X Y plane, that is why rotating with respect to Y axis
      G4double rotationAngle1 = 100*deg;
      rotm1.rotateX(coneAngle/2); // similar as the first one
      G4ThreeVector position1 = G4ThreeVector(secondMLCPos_x,secondMLCPos_y,secondMLCPos_z);
      G4Transform3D transform1 = G4Transform3D(rotm1,position1);


         secondMLC_phys = new G4PVPlacement(transform1,             //rotation,position
					    "secondMLC_phys",
	                                    secondMLCLogicV,

					    PVWorld,
					    false,
					    0,
					    false);
      rotm1.rotateX(-MLCArcAngle - coneAngle); // similar as the third one
      transform1 = G4Transform3D(rotm1,position1);
      
      fourthMLC_phys = new G4PVPlacement(transform1,             //rotation,position
				  "secondMLC_phys",
				  secondMLCLogicV,

				  PVWorld,
				  false,
				  0,
				  false);
//       
      // Region for cuts
      G4Region *regVol;
      regVol= new G4Region("MLCR");
      G4ProductionCuts* cuts = new G4ProductionCuts;
      cuts->SetProductionCut(1.0*cm);
      regVol->SetProductionCuts(cuts);
      
      firstMLCLogicV->SetRegion(regVol);
      regVol->AddRootLogicalVolume(firstMLCLogicV);
      
      secondMLCLogicV->SetRegion(regVol);
      regVol->AddRootLogicalVolume(secondMLCLogicV);
      
      
      G4VisAttributes* simpleAlSVisAtt;
      simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
      simpleAlSVisAtt->SetVisibility(true);
      firstMLCLogicV->SetVisAttributes(simpleAlSVisAtt);
      secondMLCLogicV->SetVisAttributes(simpleAlSVisAtt);


}
void  MRLinacHead::GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume)
{
          // create the accelerator-world box
      G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      G4ThreeVector halfSize;
      G4ThreeVector initialCentre;
      G4double isoCentre = 0;
      initialCentre.set(0.*mm, 0.*mm, -isoCentre);
      G4double LinacWorldDimX = 500*mm;
      G4double LinacWorldDimY = 500*mm;
      G4double LinacWorldDimZ = 1000*mm;
      halfSize.set(LinacWorldDimX,LinacWorldDimY,LinacWorldDimZ);
      G4Box *accWorldB = new G4Box("accWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
      accWorldLV = new G4LogicalVolume(accWorldB, Vacuum, "accWorldL", 0, 0, 0);
      G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
      simpleAlSVisAtt->SetVisibility(true);
// 		simpleAlSVisAtt->SetForceWireframe(false);
      accWorldLV->SetVisAttributes(simpleAlSVisAtt);
      
      PVWorld = physicalVolume;

//       PVWorld = new G4PVPlacement(0, initialCentre, "acceleratorBox", accWorldLV, PVWorld, false, 0);
      
      maxConeAngle = 8.04*2*deg;
//       diameterOfPrimaryCollimatorTop = 13.2*mm; // why it is 13.2?
       diameterOfPrimaryCollimatorTop = 7*mm;
      
      Target();
      PrimaryCollimator();
      ElectronAbsorber();
      IonizationChamber();
      PhaseSpacePlane();
//       MLC();
//       MLCSteepEdge();
      
// 	G4double fieldSize =  22.41*cm;
//       G4double fieldSize = 18.26*cm;
      G4double fieldSize = 14.11*cm;
//       G4double fieldSize = 9.96*cm;
//       G4double fieldSize = 5.82*cm;
//       G4double fieldSize = 4.1*cm;
//       G4double fieldSize = 1.66*cm;
//       
      G4double coneAngleOpen;
      coneAngleOpen = atan((fieldSize-diameterOfPrimaryCollimatorTop)/1800)/pi*180*deg;
      MLCCustomized(coneAngleOpen*2);
//       MLCCustomized(0*deg);



}




void MRLinacHead::PhaseSpacePlane()
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
  G4Transform3D transform1 = G4Transform3D(rotm1,G4ThreeVector(0,0,750*mm));
  new G4PVPlacement(transform1,"phaseSpacePlane",phaseSpacePlane_log,PVWorld,false,0,false);
  
  G4VisAttributes* simpleAlSVisAtt;
  simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
  simpleAlSVisAtt->SetVisibility(true);
  phaseSpacePlane_log->SetVisAttributes(simpleAlSVisAtt);

  


}




