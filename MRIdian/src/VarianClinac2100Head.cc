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


#include "VarianClinac2100Head.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <iostream>
using namespace std;


VarianClinac2100Head::VarianClinac2100Head()
{
}

VarianClinac2100Head::~VarianClinac2100Head()
{
}



void VarianClinac2100Head::writeInfo()
{
	std::cout <<"\n\n\tnominal beam energy: "<<idEnergy << G4endl;
	std::cout <<"\tJaw X aperture: 1) "<< jaw1XAperture/mm<<"[mm]\t2) " << jaw2XAperture/mm<< " [mm]"<< G4endl;
	std::cout <<"\tJaw Y aperture: 1) "<< jaw1YAperture/mm<<"[mm]\t2) " << jaw2YAperture/mm<< " [mm]\n"<< G4endl;
}
G4Material * VarianClinac2100Head::otherMaterials(const G4String materialName)
{
	G4Material * material=0;
	G4double A, Z, d;
	G4String name;

   // General elements
 
	A = 12.011*g/mole;
	G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);  

	A = 30.974*g/mole;
	G4Element* elP = new G4Element("Phosphorus","P",Z = 15.,A);
	 
	A = 32.064*g/mole;
	G4Element* elS = new G4Element("Sulfur","S",Z = 16.,A);
	 
	A = 55.85*g/mole;
	G4Element* elFe  = new G4Element("Iron","Fe",Z = 26.,A);
	 
	A = 51.9961*g/mole;
	G4Element* elCr = new G4Element("Cromium","Cr", Z = 24.,A);

	A = 54.94*g/mole;
	G4Element* elMn = new G4Element("Manganese","Mn", Z = 25.,A);

	A =  58.69*g/mole;
	G4Element* elNi = new G4Element("Nickel","Ni", Z = 28.,A);

	A = 28.09*g/mole;
	G4Element* elSi = new G4Element("Silicon","Si", Z = 14.,A);

	A = 183.84*g/mole;
	G4Element* elW = new G4Element("Tungsten","W", Z = 74.,A);

	if (materialName=="steel1")
	{
		d = 7.76 *g/cm3;
		G4Material* steel1 = new G4Material("steel1", d,4);
		steel1 -> AddElement(elFe, 0.935);
		steel1 -> AddElement(elS,0.01);
		steel1 -> AddElement(elMn,0.05);
		steel1 -> AddElement(elC,0.005);
		material=steel1;
	}
	else if (materialName=="steel2")
	{
		d = 8.19*g/cm3;
		G4Material* steel2 = new G4Material("steel2", d,5);
		steel2 -> AddElement(elFe, 0.759);
		steel2 -> AddElement(elNi, 0.11);
		steel2 -> AddElement(elSi,0.01);
		steel2 -> AddElement(elCr,0.12);
		steel2 -> AddElement(elP,0.001);
		material=steel2;
	}
	else if (materialName=="steel3")
	{
		d = 8.19*g/cm3;
		G4Material* steel3 = new G4Material("steel3", d,5);
		steel3 -> AddElement(elFe, 0.69);
		steel3 -> AddElement(elNi, 0.1);
		steel3 -> AddElement(elSi,0.01);
		steel3 -> AddElement(elCr,0.18);
		steel3 -> AddElement(elMn,0.02);
		material=steel3;
	}
	else if (materialName=="EZcut")
	{
		d = 7.9*g/cm3;
		G4Material* EZcut20 = new G4Material("EZcut", d,2);
		EZcut20 -> AddElement(elFe, 0.98);
		EZcut20 -> AddElement(elMn,0.02);
		material=EZcut20;
	}
	else if (materialName=="W")
	{
		d = 15*g/cm3;
		G4Material* W = new G4Material("W", d,1);
		W -> AddElement(elW, 1.);
		material=W;
	}
	return material;
}


void VarianClinac2100Head::GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume)
{
  // create the accelerator-world box
  
    PVWorld = physicalVolume;
    
    relativePositionOfLinacHead.set(0,0,450*mm);
  
    

    G4double iso = 1000; // this is the SSD distance, here we take it 100 cm
    setIsoCentre(iso);
    target();
    BeWindow();
    ionizationChamber();
    flatteningFilter();
    mirror();
    primaryCollimator();
    MLC();
    Jaw1X();
    Jaw2X();
    Jaw1Y();
    Jaw2Y();
    
    PhaseSpacePlane(200*mm);
    
    G4LogicalVolume* testLV = PVWorld->GetLogicalVolume();
    cout<<"The number of daughter logical volumes are "<<testLV->GetNoDaughters()<<endl;
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {
      cout<<"The name of the daughter logical volume is "<<testLV->GetDaughter(i)->GetName()<<endl;
    }
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)// set the position of linac components according to the needs
    {
      G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
      cout<<"The initial position information of daughter volume are as: "<<position.getX()<<" "<<position.getY()<<" "<<position.getZ()<<endl;
      position.setX(position.getX()+ relativePositionOfLinacHead.getX());
      position.setY(position.getY()+ relativePositionOfLinacHead.getY());
      position.setZ(position.getZ()+ relativePositionOfLinacHead.getZ());
      testLV->GetDaughter(i)->SetTranslation(position);
    }
    
    double maxz = testLV->GetDaughter(0)->GetTranslation().getZ();
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {
      G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
      double z = position.getZ();
      if (z>=maxz)
      {
	maxz = z;
      }
    }
    
    double minz = testLV->GetDaughter(0)->GetTranslation().getZ();
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {      
      G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
      double z = position.getZ();
      if (z<= minz )
      {
	minz = z;
      }
    }
    
    cout<<"The maxz is "<<maxz<<endl;
    cout<<"The minz is "<<minz<<endl;
    
    double zr = (maxz + minz)/2;
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {
       G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
       position.setZ(2*zr - position.getZ());
       testLV->GetDaughter(i)->SetTranslation(position);
    }
    
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {
       G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
       cout<<"The position information of daughter volume after flip are as: "<<position.getX()<<" "<<position.getY()<<" "<<position.getZ()<<endl;
       if (i == 0)
       {
	  targetPosition.setX(position.getX());
	  targetPosition.setY(position.getY());
	  targetPosition.setZ(position.getZ());
       }
    }
    
    
    for (int i = 0;i<testLV->GetNoDaughters();i++)
    {
      
      G4RotationMatrix* rotation1 = testLV->GetDaughter(i)->GetRotation();
      
      G4RotationMatrix* rotation = new G4RotationMatrix();
      rotation->rotateX(180*deg);
//       rotation->rotateY(180*deg);
//       rotation->rotateZ(90*deg);

	  
      if (rotation1)
      {
	cout<<"angle is "<<rotation1->axisAngle().getDelta()/deg<<endl;
	      cout<<"The name of the daughter logical volume is "<<testLV->GetDaughter(i)->GetName()<<endl;
	      cout<<"axis is "<<rotation1->axisAngle().axis().getX()<<" "<<rotation1->axisAngle().axis().getY()<<" "<<rotation1->axisAngle().axis().getZ()<<endl;
	      if (rotation1->axisAngle().axis().getX()!=0)
	      {
		rotation->rotateX((rotation1->axisAngle().getDelta())*rotation1->axisAngle().axis().getX()*(-1) );
	      }
	      if (rotation1->axisAngle().axis().getY()!=0)
	      {
		rotation->rotateY(rotation1->axisAngle().getDelta()*rotation1->axisAngle().axis().getY());
	      }
	      if (rotation1->axisAngle().axis().getZ()!=0)
	      {
		rotation->rotateZ(rotation1->axisAngle().getDelta()*rotation1->axisAngle().axis().getZ());
	      }
 

    }
    
    testLV->GetDaughter(i)->SetRotation(rotation);  


    }


}




void VarianClinac2100Head::reset()
{
	leavesA.clear();
	leavesB.clear();
}

void VarianClinac2100Head::target()
{
  //    materials  

    G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
    G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

	//    colors

    G4Colour  cyan    (0.0, 1.0, 1.0);
    G4Colour  magenta (1.0, 0.0, 1.0); 
 
  //    volumes
  //    beam line along z axis
  //------------------------target 6MV------------------------
    G4double targetADim_x = 0.6*cm;
    G4double targetADim_y = 0.6*cm;
    G4double targetADim_z = 0.04445*cm;
    G4Box* targetA_box = new G4Box("targetA_box",targetADim_x,targetADim_y,targetADim_z);
    G4LogicalVolume *targetA_log = new G4LogicalVolume(targetA_box,W,"targetA_log",0,0,0);
    G4double targetAPos_x = 0.0*m;
    G4double targetAPos_y = 0.0*m;
    G4double targetAPos_z = 0.20055*cm;
    targetA_phys = new G4PVPlacement(0,
	      G4ThreeVector(targetAPos_x,targetAPos_y,targetAPos_z),
	      "targetA",targetA_log,PVWorld,false,0);

    G4double targetBDim_x = 0.6*cm;
    G4double targetBDim_y = 0.6*cm;
    G4double targetBDim_z = 0.07874*cm;
    G4Box* targetB_box = new G4Box("targetB_box",targetBDim_x,targetBDim_y,targetBDim_z);
    G4LogicalVolume *targetB_log = new G4LogicalVolume(targetB_box,Cu,"targetB_log",0,0,0);
    G4double targetBPos_x = 0.0*m;
    G4double targetBPos_y = 0.0*m;
    G4double targetBPos_z = 0.07736*cm;
    targetB_phys = new G4PVPlacement(0,
	      G4ThreeVector(targetBPos_x,targetBPos_y,targetBPos_z),
	      "targetB",targetB_log,PVWorld,false,0);


	// ***********  REGIONS for CUTS

  	G4Region *regVol;
	regVol= new G4Region("targetR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	targetA_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(targetA_log);
	targetB_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(targetB_log);

	// Visualization attributes

	G4VisAttributes* simpleWSVisAtt, *simpleCuSVisAtt;
	simpleWSVisAtt= new G4VisAttributes(magenta);
	simpleWSVisAtt->SetVisibility(true);
// 	simpleWSVisAtt->SetForceSolid(true);
	simpleCuSVisAtt= new G4VisAttributes(cyan);
	simpleCuSVisAtt->SetVisibility(true);
// 	simpleCuSVisAtt->SetForceSolid(true);
	targetA_log->SetVisAttributes(simpleWSVisAtt);
	targetB_log->SetVisAttributes(simpleCuSVisAtt);

}
void VarianClinac2100Head::primaryCollimator()
{
  
  //    materials 

    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

  //    colors

   G4Colour  magenta (1.0, 0.0, 1.0); 
 
 //---------rotation matrix first collimator --------

  G4RotationMatrix*  rotateMatrix=new G4RotationMatrix();
  rotateMatrix->rotateX(180.0*deg);

  //-------------------- the first collimator upper----------------

  
  G4double innerRadiusOfTheTubeEx = 1.0*cm;
  G4double outerRadiusOfTheTubeEx = 8.*cm;
  G4double hightOfTheTubeEx = 3.0*cm;
  G4double startAngleOfTheTubeEx = 0.*deg;
  G4double spanningAngleOfTheTubeEx = 360.*deg;
  G4Tubs* UpperCollimator = new G4Tubs("UpperCollimator",innerRadiusOfTheTubeEx,
                                    outerRadiusOfTheTubeEx,hightOfTheTubeEx,
				    startAngleOfTheTubeEx,spanningAngleOfTheTubeEx);
  G4LogicalVolume *UpperCollimator_log = new G4LogicalVolume(UpperCollimator,W,"UpperCollimator_log",0,0,0);

  G4double UpperCollimatorPosX = 0.*cm;
  G4double UpperCollimatorPosY = 0.*cm;
  G4double UpperCollimatorPosZ = -1.*cm;
  UpperCollimator_phys = new G4PVPlacement(0,
					   G4ThreeVector(UpperCollimatorPosX,UpperCollimatorPosY,
							 UpperCollimatorPosZ),"UpperCollimator",
					   UpperCollimator_log,PVWorld,false,0);

 
  //-------------------- the first collimator lower----------------

  G4double  pRmin1 = 0.*cm;
 
  G4double  pRmax1 = 0.5*cm;
  G4double  pRmin2 = 0.*cm;
  G4double  pRmax2 = 1.7658592*cm;
  G4double  hightOfTheCone =3.2*cm;
  G4double  startAngleOfTheCone = 0.*deg;
  G4double  spanningAngleOfTheCone = 360.*deg;

  G4Cons* collim_cone = new G4Cons("collim_cone",pRmin1,pRmax1,pRmin2,
				   pRmax2,hightOfTheCone,startAngleOfTheCone,
				   spanningAngleOfTheCone);
  G4LogicalVolume *collim_log = new G4LogicalVolume(collim_cone,Vacuum,"collim_log",0,0,0);


  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 8.*cm;
  G4double hightOfTheTube = 3.1*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
                                    outerRadiusOfTheTube,hightOfTheTube,
				    startAngleOfTheTube,spanningAngleOfTheTube);
//  G4LogicalVolume *tracker_log = new G4LogicalVolume(tracker_tube,W,"tracker_log",0,0,0);


  G4SubtractionSolid* CylMinusCone = new G4SubtractionSolid("Cyl-Cone",
  							tracker_tube,collim_cone);
  G4LogicalVolume *CylMinusCone_log = new G4LogicalVolume(CylMinusCone,W,"CylminusCone_log",0,0,0);
  G4double CminusCPos_x = 0.*cm;
  G4double CminusCPos_y = 0.*cm;
  G4double CminusCPos_z = +6.2*cm;
  CylMinusCone_phys = new G4PVPlacement(rotateMatrix,
					G4ThreeVector(CminusCPos_x,CminusCPos_y,CminusCPos_z),
					"CylMinusCone",CylMinusCone_log,PVWorld,false,0);
  
//--------- Visualization attributes -------------------------------
   G4VisAttributes* simpleTungstenWVisAtt= new G4VisAttributes(magenta);
   simpleTungstenWVisAtt->SetVisibility(true);
//    simpleTungstenWVisAtt->SetForceSolid(true);
   collim_log->SetVisAttributes(simpleTungstenWVisAtt);
   
 
   CylMinusCone_log->SetVisAttributes(simpleTungstenWVisAtt);
   UpperCollimator_log->SetVisAttributes(simpleTungstenWVisAtt);

	// ***********  REGIONS for CUTS

   G4Region *regVol;
	regVol= new G4Region("PrymCollR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	collim_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(collim_log);

	CylMinusCone_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(CylMinusCone_log);

	UpperCollimator_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(UpperCollimator_log);

}
void VarianClinac2100Head::BeWindow()
{
	G4Material *Be=G4NistManager::Instance()->FindOrBuildMaterial("G4_Be");
	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAtt;
		// Region for cuts
	regVol= new G4Region("BeWindow");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	G4Tubs* BeWTube = new G4Tubs("BeWindowTube", 0., 36.*mm, 0.2*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *BeWTubeLV = new G4LogicalVolume(BeWTube, Be, "BeWTubeLV", 0, 0, 0);
        BeWTubePV=new G4PVPlacement(0, G4ThreeVector(0.,0.,100.*mm), "BeWTubePV", BeWTubeLV,
                                    PVWorld, false, 0);

	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	BeWTubeLV->SetVisAttributes(simpleAlSVisAtt);
	BeWTubeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(BeWTubeLV);

}
void VarianClinac2100Head::flatteningFilter()
{

	G4double z0, h0;
	G4ThreeVector centre, halSize;
	G4Material *Cu=G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region("flatfilterR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.5*cm);
	regVol->SetProductionCuts(cuts);

	G4VisAttributes* simpleAlSVisAtt;

	// one
	z0=130.0*mm;
	h0=5.0/2.*cm;
	centre.set(0.,0.,z0);
	G4Cons *FFL1A_1Cone = new G4Cons("FFL1A_1", 0.*cm, 0.3*cm, 0.*cm, 5.*cm, h0, 0.*deg, 360.*deg);
	G4LogicalVolume *FFL1A_1LV = new G4LogicalVolume(FFL1A_1Cone, Cu, "FFL1A_1LV", 0, 0, 0);
	FFL1A_1PV=new G4PVPlacement(0, centre, "FFL1A_1PV", FFL1A_1LV, PVWorld, false, 0);

	// two
	z0+=h0;
	h0=0.081/2.*cm;
	z0+=h0;
	centre.setZ(z0);
	z0+=h0;
	G4Tubs *FFL2_1Tube = new G4Tubs("FFL6_1", 0.*cm, 2.5*cm, h0, 0.*deg, 360.*deg);
	G4LogicalVolume *FFL2_1LV = new G4LogicalVolume(FFL2_1Tube, Cu, "FFL2_1LV", 0, 0, 0);
	FFL2_1PV=new G4PVPlacement(0, centre, "FFL2_1PV", FFL2_1LV, PVWorld, false, 0);

	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	FFL1A_1LV->SetVisAttributes(simpleAlSVisAtt);
	FFL2_1LV->SetVisAttributes(simpleAlSVisAtt);

	FFL1A_1LV->SetRegion(regVol);
	FFL2_1LV->SetRegion(regVol);

	regVol->AddRootLogicalVolume(FFL1A_1LV);
	regVol->AddRootLogicalVolume(FFL2_1LV);
}
void VarianClinac2100Head::ionizationChamber()
{


	G4Material *material=G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	G4VisAttributes* simpleAlSVisAtt;
		// Region for cuts
	G4Region *regVol;
	regVol= new G4Region("ionizationChamber");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	G4Tubs* ICTubeW = new G4Tubs("ionizationChamberTube", 0., 2.*2.54*10.*mm, 0.016*25.4*mm, 0.*deg, 360.*deg);
	G4Tubs* ICTubeP = new G4Tubs("ionizationChamberTube", 0., 2.*2.54*10.*mm, 0.010*25.4*mm, 0.*deg, 360.*deg);

	G4ThreeVector centre;
	// W1
	centre.set(0.,0.,157.*mm);
	G4LogicalVolume *PCUTubeW1LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW1LV", 0, 0, 0);
	PCUtubeW1PV=new G4PVPlacement(0, centre, "ionizationChamberTubeW1PV", PCUTubeW1LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeW1LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW1LV);

	// P1
	centre.set(0.,0.,158.*mm);
	G4LogicalVolume *PCUTubeP1LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP1LV", 0, 0, 0);
	PCUtubeP1PV=new G4PVPlacement(0, centre, "ionizationChamberTubeP1PV", PCUTubeP1LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeP1LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP1LV);

	// W2
	centre.set(0.,0.,159.*mm);
	G4LogicalVolume *PCUTubeW2LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW2LV", 0, 0, 0);
	PCUtubeW2PV=new G4PVPlacement(0, centre, "ionizationChamberTubeW2PV", PCUTubeW2LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeW2LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW2LV);

	// P2
	centre.set(0.,0.,160.*mm);
	G4LogicalVolume *PCUTubeP2LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP2LV", 0, 0, 0);
	PCUtubeP2PV=new G4PVPlacement(0, centre, "ionizationChamberTubeP2PV", PCUTubeP2LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeP2LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP2LV);

	// W3
	centre.set(0.,0.,161.*mm);
	G4LogicalVolume *PCUTubeW3LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW3LV", 0, 0, 0);
	PCUtubeW3PV=new G4PVPlacement(0, centre, "ionizationChamberTubeW3PV", PCUTubeW3LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeW3LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW3LV);

	// P3
	centre.set(0.,0.,162.*mm);
	G4LogicalVolume *PCUTubeP3LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP3LV", 0, 0, 0);
	PCUtubeP3PV=new G4PVPlacement(0, centre, "ionizationChamberTubeP3PV", PCUTubeP3LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	PCUTubeP3LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP3LV);

}
void VarianClinac2100Head::mirror()
{
	G4Material *MYLAR=G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
	G4VisAttributes* simpleAlSVisAtt;
		// Region for cuts
	G4Region *regVol;
	regVol= new G4Region("Mirror");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	G4Tubs* MirrorTube = new G4Tubs("MirrorTube", 0., 63.*mm, .5*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *MirrorTubeLV = new G4LogicalVolume(MirrorTube, MYLAR, "MirrorTubeLV", 0, 0, 0);
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	cRotation->rotateY(12.0*deg);
        MirrorTubePV=new G4PVPlacement(cRotation, G4ThreeVector(0., 0., 175.*mm), "MirrorTubePV", MirrorTubeLV,PVWorld, false, 0);

	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	MirrorTubeLV->SetVisAttributes(simpleAlSVisAtt);
	MirrorTubeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(MirrorTubeLV);

}

void VarianClinac2100Head::SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4double aperture, G4RotationMatrix *cRotation)
{
  
    cout<<"idJaw = "<<idJaw<<" halfSize = "<<halfSize<<" aperture = "<<aperture<<endl;
	using namespace std;
        G4double theta, x, y, z, dx, dy, dz;
	x=centre.getX();
	y=centre.getY();
	z=centre.getZ();
	theta=fabs(atan(aperture/isoCentre));
	
	cout<<"theta in deg = "<<theta/deg<<endl;
	
	dx=halfSize.getX();
	dy=halfSize.getY();
	dz=halfSize.getZ();
	
//        G4double p1x, p1y, p2x, p2y;
	
	switch (idJaw)
	{
	case 1: //idJaw1XV2100:
		centre.set(z*sin(theta)+dx*cos(theta), y, z*cos(theta)-dx*sin(theta));
//		beta=fabs(atan(dx/dz));
//		R=std::sqrt(dx*dx+dz*dz);
//		p1x=centre.getX()-R*sin(theta+beta);
//		p1y=centre.getZ()-R*cos(theta+beta);
//		p2x=centre.getX()+R*sin(theta-beta);
//		p2y=centre.getZ()+R*cos(theta-beta);
		
		cRotation->rotateY(-theta);
		halfSize.set(fabs(dx*cos(theta)+dz*sin(theta)), fabs(dy), fabs(dz*cos(theta)+dx*sin(theta)));
		break;
	case 2: //idJaw2XV2100:
		centre.set(-(z*sin(theta)+dx*cos(theta)), y, z*cos(theta)-dx*sin(theta));
//		beta=fabs(atan(dx/dz));
//		R=std::sqrt(dx*dx+dz*dz);
//		p1x=centre.getX()+R*sin(theta+beta);
//		p1y=centre.getZ()-R*cos(theta+beta);
//		p2x=centre.getX()-R*sin(theta-beta);
//		p2y=centre.getZ()+R*cos(theta-beta);
		
		cRotation->rotateY(theta);
		halfSize.set(fabs(dx*cos(theta)+dz*sin(theta)), fabs(dy), fabs(dz*cos(theta)+dx*sin(theta)));
		break;
	case 3: //idJaw1YV2100:
		centre.set(x, z*sin(theta)+dy*cos(theta), z*cos(theta)-dy*sin(theta));
//		beta=fabs(atan(dy/dz));
//		R=std::sqrt(dy*dy+dz*dz);
//		p1x=centre.getY()-R*sin(theta+beta);
//		p1y=centre.getZ()-R*cos(theta+beta);
//		p2x=centre.getY()+R*sin(theta-beta);
//		p2y=centre.getZ()+R*cos(theta-beta);
		
		cRotation->rotateX(theta);
		halfSize.set(fabs(dx), fabs(dy*cos(theta)+dz*sin(theta)), fabs(dz*cos(theta)+dy*sin(theta)));
		break;
	case 4: //idJaw2YV2100:
		centre.set(x, -(z*sin(theta)+dy*cos(theta)), z*cos(theta)-dy*sin(theta));
//		beta=fabs(atan(dy/dz));
//		R=std::sqrt(dy*dy+dz*dz);
//		p1x=centre.getY()+R*sin(theta+beta);
//		p1y=centre.getZ()-R*cos(theta+beta);
//		p2x=centre.getY()-R*sin(theta-beta);
//		p2y=centre.getZ()+R*cos(theta-beta);
		
		cRotation->rotateX(-theta);
		halfSize.set(fabs(dx), fabs(dy*cos(theta)+dz*sin(theta)), fabs(dz*cos(theta)+dy*sin(theta)));
		break;
	}
}



void VarianClinac2100Head::Jaw1X()
{
	G4Material *steel1=otherMaterials("steel1");
	G4String name="Jaws1X";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(320.+80./2.)*mm);
	halfSize.set(45.*mm, 93.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	
	setJaw1X(-100*mm);
	
	SetJawAperture(1, centre, halfSize, jaw1XAperture, cRotation);
	
	cout<<"cRotation information "<<cRotation->getTheta()<<" "<<cRotation->getPhi()<<" "<<cRotation->getPsi()<<endl;
        phVol1X= new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

}
void VarianClinac2100Head::Jaw2X()
{
	G4Material *steel1=otherMaterials("steel1");
	G4String name="Jaws2X";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(320.+80./2.)*mm);
	halfSize.set(45.*mm, 93.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	
	setJaw2X(100*mm);
	SetJawAperture(2, centre, halfSize, jaw2XAperture, cRotation);
        phVol2X= new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

}
void VarianClinac2100Head::Jaw1Y()
{
	G4Material *steel1=otherMaterials("steel1");
	G4String name="Jaws1Y";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(230.+80./2.)*mm);
	halfSize.set(93.*mm, 35.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	
	setJaw1Y(-100*mm);
	
	SetJawAperture(3, centre, halfSize, jaw1YAperture, cRotation);
        phVol1Y= new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

}
void VarianClinac2100Head::Jaw2Y()
{
	G4Material *steel1=otherMaterials("steel1");
	G4String name="Jaws2Y";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(230.+80./2.)*mm);
	halfSize.set(93.*mm, 35.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	
	setJaw2Y(100*mm);
	
	SetJawAperture(4, centre, halfSize, jaw2YAperture, cRotation);
        phVol2Y= new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts

	G4Region *regVol;
	regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

}
void VarianClinac2100Head::MLC()
{

 //    material
	G4Material *Fe=G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region("MLCR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.0*cm);
	regVol->SetProductionCuts(cuts);

	G4ThreeVector boxSize;

	G4ThreeVector centreStart;
	centreStart.set(0.,0.,(330.+600.)/2.*mm);

	boxSize.set(6./2.*mm, 180./2.*mm, 50./2.*mm);

	// single leaf
	G4Box* boxLeaf =new G4Box("LeafBox", boxSize.getX(), boxSize.getY(), boxSize.getZ());

	G4LogicalVolume *leafLVA = new G4LogicalVolume(boxLeaf, Fe, "leafSolidALV", 0, 0, 0);
	G4LogicalVolume *leafLVB = new G4LogicalVolume(boxLeaf, Fe, "leafSolidBLV", 0, 0, 0);

	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	leafLVA->SetVisAttributes(simpleAlSVisAtt);
	leafLVA->SetRegion(regVol);
	regVol->AddRootLogicalVolume(leafLVA);
	
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
// 	simpleAlSVisAtt->SetForceSolid(true);
	leafLVB->SetVisAttributes(simpleAlSVisAtt);
	leafLVB->SetRegion(regVol);
	regVol->AddRootLogicalVolume(leafLVB);

	int i;
	G4String PVname;
	int j=0;

	G4ThreeVector centre;
	
	
	for(int i = 0;i<10; i++)
	{
	  setLeavesAx(200*mm);
	  setLeavesBx(200*mm);
	}
	
	
	for (int i = 0;i<leavesA.size();i++)
	{
	  cout<<"leavesA "<<i<<" "<<leavesA[i]<<endl;
	}
	
	for (int i = 0;i<leavesB.size();i++)
	{
	  cout<<"leavesA "<<i<<" "<<leavesB[i]<<endl;
	}
	

	int nhalfLeaves=(int)(leavesA.size()/2.);
	centre= centreStart + G4ThreeVector(-nhalfLeaves*boxSize.getX(), 0.,0.);
	
		

	for (i=1;i<(int)leavesA.size(); i++)
	{
		G4String str;
		char appo[12];
		sprintf(appo,"%d",i);
		str=appo;
		PVname="leafA"+str;
		centre.setX(centre.getX()+boxSize.getX()*2.);
		centre.setY(-boxSize.getY()-leavesA[i]);
		leafPhys=new G4PVPlacement(0, centre, PVname, leafLVA, PVWorld, false, i);
		j++;
	}
	nhalfLeaves=(int)(leavesB.size()/2.);
	centre=centreStart+G4ThreeVector(-nhalfLeaves*boxSize.getX(), 0.,0.);
	for (i=1;i<(int)leavesB.size(); i++)
	{
		G4String str;
		char appo[12];
		sprintf(appo,"%d",i);
		str=appo;
		PVname="leafB"+str;
		centre.setX(centre.getX()+boxSize.getX()*2.);
		centre.setY(+boxSize.getY()+leavesB[i]);
		leafPhys=new G4PVPlacement(0, centre, PVname, leafLVB, PVWorld, false, i);
		j++;
	}
}


void VarianClinac2100Head::PhaseSpacePlane(G4double zPosition)
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
  G4Transform3D transform1 = G4Transform3D(rotm1,G4ThreeVector(0,0,zPosition));
  new G4PVPlacement(transform1,"phaseSpacePlane",phaseSpacePlane_log,PVWorld,false,0,false);
  
  G4VisAttributes* simpleAlSVisAtt;
  simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
  simpleAlSVisAtt->SetVisibility(true);
  phaseSpacePlane_log->SetVisAttributes(simpleAlSVisAtt);

}

