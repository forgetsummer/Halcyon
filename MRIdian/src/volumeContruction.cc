#include "volumeConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4TransportationManager.hh"
#include "G4RegionStore.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Region.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Paraboloid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

volumeConstruction::volumeConstruction()

: pMere(0),vacum(0),
WinTitanium(0),Target(0),CollimatorA1(0),CollimatorA2(0),ColliSec(0),LoGiCPlqAl(0),LoGiCJawX1_1(0),
LoGiCJawX1_2(0),LoGiCJawX1_3(0),LoGiCJawX1_4(0),LoGiCJawX1_5(0),LoGiCJawX1_6(0),LoGiCJawY1_1(0),LoGiCJawY1_2(0),
LoGiCJawY1_3(0),LoGiCJawY1_4(0),LoGiCJawY1_5(0),pvacum(0),pWinTitanium(0)
,pTarget(0),pCollimatorA1(0),pCollimatorA2(0),pColliSec(0),pLoGiCPlqAl(0), pLoGiCJawX1_1(0),pLoGiCJawX1_2(0),
pLoGiCJawX1_3(0),pLoGiCJawX1_4(0),pLoGiCJawX1_5(0),pLoGiCJawX1_6(0),pLoGiCJawY1_1(0),pLoGiCJawY1_2(0),
pLoGiCJawY1_3(0),pLoGiCJawY1_4(0),pLoGiCJawY1_5(0),pLoGiCJawX2_1(0),pLoGiCJawX2_2(0),
pLoGiCJawX2_3(0),pLoGiCJawX2_4(0),pLoGiCJawX2_5(0),pLoGiCJawX2_6(0),pLoGiCJawY2_1(0),pLoGiCJawY2_2(0),
pLoGiCJawY2_3(0),pLoGiCJawY2_4(0),pLoGiCJawY2_5(0),
pChI_ElementTypeB2_layer1(0),pChI_ElementTypeA2_layer2(0),pChI_ElementTypeB2_layer3(0),
pChI_ElementTypeB2_layer4(0),pChI_ElementTypeA2_layer5(0),pChI_ElementTypeB2_layer6(0)


{


dec=-25*cm;


  }
// ==================================================================================
//  DESCRECTOR METHOD OF THE CLASS VOLUMECONSTRUCTION
// ==================================================================================

volumeConstruction::~volumeConstruction()
{; }
// ==================================================================================
//  THIS METHOD SET-UP VISUALIZATION ATTRIBUT FOR EACH ELEMENT OF ACCELERATOR
// ==================================================================================

 G4VisAttributes *volumeConstruction::AttributName(const G4String attributName)
{
G4VisAttributes * attribut=new G4VisAttributes();
if (attributName=="Titanium_attribut")
	{
  G4VisAttributes * Titanium_attribut=  new G4VisAttributes( G4Colour(255/255. ,0/255. ,255/255. )); 
 attribut=Titanium_attribut;
return attribut;
}
if (attributName=="phantom_attribut")
	{
  G4VisAttributes * phantom_attribut=  new G4VisAttributes( G4Colour(127/255. ,0/255. ,255/255. )); 
 attribut=phantom_attribut;
return attribut;
}
if (attributName=="water_attribut")
	{
  G4VisAttributes * water_attribut=  new G4VisAttributes( G4Colour(0/255. ,0/255. ,255/255. )); 
 attribut=water_attribut;
return attribut;
}
if (attributName=="Tungsten_attribut")
	{
  G4VisAttributes * Tungsten_attribut=   new G4VisAttributes( G4Colour(15/255. ,155/255. ,0/255. ));
 attribut=Tungsten_attribut;
return attribut;
}

if (attributName=="Pb_attribut")
	{
  G4VisAttributes * Pb_attribut=   new G4VisAttributes( G4Colour(43/255. ,21/255. ,145/255. ));
attribut=Pb_attribut;
return attribut;
}

if (attributName=="WNICU_attribut")
	{
  G4VisAttributes * WNICU_attribut=  new   G4VisAttributes( G4Colour(255/255. ,21/255. ,0/255. ));
attribut=WNICU_attribut;
return attribut;
}

if (attributName=="Stainless_steel_attribut")
	{
  G4VisAttributes * Stainless_steel_attribut=   new G4VisAttributes( G4Colour(204/255. ,0/255. ,204/255. ));
attribut=Stainless_steel_attribut;
return attribut;
}


if (attributName=="XC10_attribut")
	{
  G4VisAttributes * XC10_attribut=   new G4VisAttributes( G4Colour(0/255. ,255/255. ,0/255. ));
attribut= XC10_attribut;
return attribut;

}
if (attributName=="Al_attribut")
	{
  G4VisAttributes * Al_attribut= new G4VisAttributes( G4Colour(255/255. ,255/255. ,0/255. ));
attribut=  Al_attribut;
return attribut;
}

if (attributName=="Kapton_attribut")
	{
  G4VisAttributes * Kapton_attribut= new G4VisAttributes( G4Colour(255/255. ,255/255. ,0/255. ));
attribut=  Kapton_attribut;
return attribut;
}
return attribut;
}

// ==================================================================================
//  THIS METHOD SET-UP MATERIAL FOR EACH ELEMENT OF ACCELERATOR
// ==================================================================================


G4Material * volumeConstruction::MaterialName(const G4String materialName)
{
G4double z, a, density;
 G4String name, symbol;
 G4int ncomponents;
G4double fractionmass;
G4Material * material=0;

	
	//--------- les Elements . -----------------------------------------------

  G4Element* elementH = new G4Element( "Hydrogen", "H", 1. , 1.00794*g/mole );
  G4Element* elementC = new G4Element( "Carbon", "C", 6. , 12.01*g/mole );
  G4Element* elementN = new G4Element( "Nitrogen", "N", 7. , 14.00674*g/mole );
  G4Element* elementO = new G4Element( "Oxygen", "O", 8. , 15.9994*g/mole );
  //G4Element* elementMg = new G4Element( "Magnesium", "Mg", 12. , 24.305*g/mole );
  G4Element* elementSi = new G4Element( "Silicon", "Si", 14. , 28.0855*g/mole );
  //G4Element* elementP = new G4Element( "Phosphorus", "P", 15. , 30.973762*g/mole );
 // G4Element* elementCl = new G4Element( "Chlorine", "Cl", 17. , 35.4527*g/mole );
  G4Element* elementAr = new G4Element( "Argon", "Ar", 18. , 39.948*g/mole );
  //G4Element* elementCa = new G4Element( "Calcium", "Ca", 20. , 40.078*g/mole );
  G4Element* elementCr = new G4Element( "Chromium", "Cr", 24. , 51.9961*g/mole );
  G4Element* elementMn = new G4Element( "Maganese", "Mn", 25. , 54.93805*g/mole );
  G4Element* elementFe = new G4Element( "Fer", "Fe", 26. , 55.85*g/mole );
  G4Element* elementNi = new G4Element( "Nickel", "Ni", 28. , 58.69*g/mole );
  G4Element* elementCu = new G4Element( "Copper", "Cu", 29. , 63.546*g/mole );
  G4Element* elementW = new G4Element( "W", "W", 74. , 183.85*g/mole );

	if (materialName=="air")
	{
  G4Material* air = new G4Material("air",density=0.0012*g/cm3,  ncomponents=4);
  air->AddElement( elementC, 0.0124 *perCent);
  air->AddElement( elementN, 75.527 *perCent);
  air->AddElement( elementO, 23.178 *perCent);
  air->AddElement( elementAr, 1.2827 *perCent);
  material=air;	
}
	else if (materialName=="Al")
	{
		G4Material* Al = new G4Material("Aliminuim", z=13., a=26.981539*g/mole, density= 2.7*g/cm3);
		material=Al;
	}
	else if (materialName=="Pb")
	{
		G4Material* Pb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.34*g/cm3);
		material=Pb;
	}
else if (materialName=="water")
	{
		G4Material* water = new G4Material("water",  density=0.9982*g/cm3,ncomponents= 2);
  water->AddElement( elementH, 11.11*perCent);
  water->AddElement( elementO, 88.89*perCent );
		 material=water;
	}
else if (materialName=="PMMA")
	{
		G4Material* PMMA = new G4Material(name="PMMA",density=1.19*g/cm3, ncomponents=3);
  PMMA->AddElement(elementH,8*perCent);
  PMMA->AddElement(elementC,60*perCent);
  PMMA->AddElement(elementO,32*perCent);
  
		 material=PMMA;
	}
	else if (materialName=="Tungsten")
	{
		 G4Material* Tungsten = new G4Material("Tungsten", z=74, a=183.85*g/mole, density=19.3*g/cm3);
		 material=Tungsten;
	}
	else if (materialName=="Titanium")
	{
		 G4Material* Titanium = new G4Material("Titanium", z=22,  a=28*g/mole, density=4.54*g/cm3 );

		material=Titanium;
	}

else if (materialName=="XC10")
	{
		G4Material* XC10 = new G4Material(name="XC10",density=7.8*g/cm3, ncomponents=3);
  XC10->AddElement(elementC, fractionmass=0.1*perCent);
  XC10->AddElement(elementMn, fractionmass=0.6*perCent);
  XC10->AddElement(elementFe, fractionmass=99.3*perCent);
		material=XC10;
	}

else if (materialName=="WNICU")
	{
		G4Material* WNICU = new G4Material(name="WNICU",density=16.8*g/cm3, ncomponents=3);
  WNICU->AddElement(elementNi, fractionmass=7*perCent);
  WNICU->AddElement(elementCu, fractionmass=2.5*perCent);
  WNICU->AddElement(elementW, fractionmass=90.5*perCent);
		material=WNICU;
	}
else if (materialName=="Stainless_steel")
	{
		G4Material* Stainless_steel = new G4Material(name="Stainless_steel",density=7.8*g/cm3, ncomponents=6);
  Stainless_steel->AddElement(elementC,0.02*perCent);
  Stainless_steel->AddElement(elementMn,2*perCent);
  Stainless_steel->AddElement(elementFe,68.98*perCent);
  Stainless_steel->AddElement(elementSi,1*perCent);
  Stainless_steel->AddElement(elementCr,18*perCent);
  Stainless_steel->AddElement(elementNi,10*perCent);
		material=Stainless_steel;
	}

else if (materialName=="Kapton")
	{
		G4Material* Kapton = new G4Material(name="Kapton",density=1.42*g/cm3, ncomponents=4);
  Kapton->AddElement(elementH,2.6360*perCent);
  Kapton->AddElement(elementC,69.113*perCent);
  Kapton->AddElement(elementN,7.327*perCent);
  Kapton->AddElement(elementO,20.924*perCent);

		material=Kapton;
	}

	return material;
}
G4VPhysicalVolume* volumeConstruction::Construct( )
{

// ==================================================================================
//  THIS METHOD SET-UP THE WORLD VOLUME
// ==================================================================================


G4Material *air=this->MaterialName("air");
//G4VisAttributes *Invisible_attribut=this->AttributName("Invisible_attribut");

   G4Box *solidvolMere= new G4Box("solidvolMere",20.05*cm, 20.05*cm, 0.3*m );
  G4LogicalVolume * volMere = new G4LogicalVolume(solidvolMere, 	 //its solid
			 air, 		 //its material
			"volMere" ,		 //its name
			 0,0,0);
//cet volume est invisible
 volMere ->SetVisAttributes(G4VisAttributes::GetInvisible());
   
G4RotationMatrix rotMatrixthispMere;   // unit rotation matrix
G4double anglethispMere = 0.0*deg;   // rotational angle
rotMatrixthispMere.rotateX(anglethispMere);  // rot matrix

  this->pMere= new G4PVPlacement(G4Transform3D(rotMatrixthispMere,	//rotation 
		 G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm)),
		 "pMere",   //its name  (2nd constructor)
		 volMere,         //its logical volume 
		 NULL,              //its mother volume 
		 false,                 //no boolean operation 
	
	 0);                       //


        WinTitaniumConstructor();
	targetConstructor();
	primaryCollimatorConstructor();
	secondaryCollimatorConstructor();
	flatteningFilteConstructor();
        ionizationChamberConstructor();
	AlPlaqueConstructor();
	jawsConstructor();
        

 return this->pMere; 
        

}



// ==================================================================================
//  THIS METHOD SET-UP THE WINDOW TITANIUM
// =================================================================================

void volumeConstruction::WinTitaniumConstructor()
{
G4VisAttributes *Titanium_attribut=this->AttributName("Titanium_attribut");
 G4Material *Titanium=this->MaterialName("Titanium");	 

//-------------	SOLID VOLUME OF WINDOW  ----------------------------
G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;

G4double Rayon_WinTetanium =4.*mm;
G4double DemiZ_WinTetanium =0.025*mm;

 G4Tubs *solidWinTetanium= new G4Tubs("solidWinTetanium", 0.0*cm, Rayon_WinTetanium, DemiZ_WinTetanium, Angle_Debut, Angle_Fin );
//-------------	LOGICAL VOLUME OF WINDOW  ----------------------------
 WinTitanium = new G4LogicalVolume(solidWinTetanium, 	 //its solid
			 Titanium, 		 //its material
			"WinTitanium" ,		 //its name
			 0,0,0);

WinTitanium->SetVisAttributes(Titanium_attribut);

//-------------	PHYSICAL VOLUME OF WINDOW  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepWinTitanium = 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepWinTitanium);  // rot matrix
G4double WinTitanium_centre_pos_z= dec-15*mm;
pWinTitanium= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, WinTitanium_centre_pos_z)),
		 "pWinTitanium",   //its name  (2nd constructor)
		 WinTitanium,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);                       //copy number 

G4Region *regWinTitanium;
       regWinTitanium= new G4Region("WinTitaniumRegion");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(0.339*mm,"gamma");
        cuts1->SetProductionCut(0.067*mm,"e-");
        cuts1->SetProductionCut(0.067*mm,"e+");
        regWinTitanium->SetProductionCuts(cuts1);

      WinTitanium->SetRegion(regWinTitanium);
        regWinTitanium->AddRootLogicalVolume(WinTitanium);

    
//-------------------------------------------------------

} 
// ==================================================================================
//  THIS METHOD SET-UP THE TARGET
// ==================================================================================
void volumeConstruction::targetConstructor()
{

//-------------	ATTRIBUT AND MATERIAL OF TARGET  ----------------------------
 G4VisAttributes *Tungsten_attribut=this->AttributName("Tungsten_attribut");
 G4Material *Tungsten=this->MaterialName("Tungsten");

//-------------	SOLID VOLUME OF TARGET  ----------------------------
//targetA
  G4double TargetA_DemiX =5.0*mm;
  G4double TargetA_DemiY =5.0*mm;
  G4double TargetA_DemiZ =7.5*mm ;
  G4double Angle_Debut =0.0*rad;
  G4double Angle_Fin =360.0*rad;

  G4Box *solidTargetA= new G4Box("solidTargetA",TargetA_DemiX,TargetA_DemiY,TargetA_DemiZ);

//targetB
  G4double TargetB_Rayon =3.0*mm;
  G4double TargetB_DemiZ =5.5*mm;
  G4Tubs *solidTargetB= new G4Tubs("solidTargetB",0.0*mm, TargetB_Rayon,TargetB_DemiZ,Angle_Debut,Angle_Fin );

//targetB et targetA
  G4ThreeVector  Target_B_to_A_decalage_z(0,0,2.0*mm);
  G4RotationMatrix* Target_B_to_A_rotation_z = new G4RotationMatrix;
  Target_B_to_A_rotation_z->rotateZ(0.*rad);

  G4SubtractionSolid* solidTarget = new G4SubtractionSolid("solidTarget",solidTargetA,solidTargetB,Target_B_to_A_rotation_z,Target_B_to_A_decalage_z);

//-------------	LOGICAL VOLUME OF TARGET  ----------------------------
  Target = new G4LogicalVolume(solidTarget, 	 //its solid
			  Tungsten, 		 //its material
			"Target" ,		 //its name
			 0,0,0);


  Target->SetVisAttributes(Tungsten_attribut);

//-------------	PHYSICAL VOLUME OF TARGET  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix
G4double Target_centre_pos_z=dec-3.5*mm;
  pTarget= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m,Target_centre_pos_z)),
		 "pTarget",   //its name  (2nd constructor)
	Target,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);  

G4Region *regTarget;
       regTarget= new G4Region("TargetRegion");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(0.0148*mm,"gamma");
        cuts1->SetProductionCut(0.037*mm,"e-");
        cuts1->SetProductionCut(0.037*mm,"e+");
        regTarget->SetProductionCuts(cuts1);

       Target->SetRegion(regTarget);
        regTarget->AddRootLogicalVolume(Target);

    
}
// ==================================================================================
//  THIS METHOD SET-UP THE PRIMARY COLLIMATOR
// ==================================================================================
void volumeConstruction::primaryCollimatorConstructor()
{


//-------------	ATTRIBUT AND MATERIAL OF PRIMARY COLLIMATOR  ----------------------------
 G4VisAttributes *XC10_attribut=this->AttributName("XC10_attribut");
 G4Material *XC10=this->MaterialName("XC10");
 G4VisAttributes *WNICU_attribut=this->AttributName("WNICU_attribut");
 G4Material *WNICU=this->MaterialName("WNICU");

//Collimateur primaire

//collimatorA1
//-------------	SOLID VOLUME OF COLLIMATOR_A1  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix

G4double ColliA1B1_Rayon =17.0*mm;
G4double ColliA1B1_DemiZ =7.93*mm;
G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;

  G4Tubs *solidColliA1B1= new G4Tubs("solidColliA1B1",0.0*mm , ColliA1B1_Rayon,ColliA1B1_DemiZ ,Angle_Debut, Angle_Fin  );

  G4double ColliA1B2_Rayon1 =17.0*mm;
  G4double ColliA1B2_Rayon2 =27.0*mm;
  G4double ColliA1B2_DemiZ =31.57*mm;

  G4Paraboloid *solidColliA1B2= new G4Paraboloid("solidColliA1B2",ColliA1B2_DemiZ,ColliA1B2_Rayon1,ColliA1B2_Rayon2);

  G4double ColliA1_Rayon =70.50*mm;
  G4double ColliA1_DemiZ =39.5*mm;
  G4Tubs *solidColliA1= new G4Tubs("solidColliA1",0.0*mm ,ColliA1_Rayon, ColliA1_DemiZ,Angle_Debut, Angle_Fin  );

  G4ThreeVector  ColliA1_B1_to_B2_decalage_z(0,0,-ColliA1_DemiZ);
  G4RotationMatrix* ColliA1_B1_to_B2_Rotation_z = new G4RotationMatrix;
  ColliA1_B1_to_B2_Rotation_z->rotateX(0.*deg);

  G4UnionSolid* CollimatorA1B1 = new G4UnionSolid("solidColliA1B",solidColliA1B2,solidColliA1B1,ColliA1_B1_to_B2_Rotation_z,ColliA1_B1_to_B2_decalage_z);


  G4ThreeVector  ColliA1_A1_to_A1B1_decalage_z(0,0,-7.93*mm);
  G4RotationMatrix* ColliA1_A1_to_A1B1_Rotation_z = new G4RotationMatrix;
  ColliA1_A1_to_A1B1_Rotation_z->rotateX(180.*deg);
  G4SubtractionSolid* SolidCollimatorA1 = new G4SubtractionSolid("SolidCollimatorA1",solidColliA1,CollimatorA1B1,ColliA1_A1_to_A1B1_Rotation_z,ColliA1_A1_to_A1B1_decalage_z);

//-------------	LOGICAL VOLUME OF COLLIMATOR_A1  ----------------------------
  CollimatorA1 = new G4LogicalVolume(SolidCollimatorA1, 	 //its solid
			  XC10, 		 //its material
			"CollimatorA1" ,		 //its name
			 0,0,0);



  CollimatorA1->SetVisAttributes(XC10_attribut);
G4Region *regCollimatorA1;
   regCollimatorA1= new G4Region("CollimatorA1Region");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(0.25*mm,"gamma");
        cuts1->SetProductionCut(0.068*mm,"e-");
        cuts1->SetProductionCut(0.068*mm,"e+");
        regCollimatorA1->SetProductionCuts(cuts1);

     CollimatorA1->SetRegion(regCollimatorA1);
        regCollimatorA1->AddRootLogicalVolume(CollimatorA1);

//collimatorA2
//-------------	SOLID VOLUME OF COLLIMATOR_A2  ----------------------------
G4double ColliA2B2_Rayon1 =26.5*mm;
G4double ColliA2B2_Rayon2 =40.5*mm;
G4double ColliA2B2_DemiZ =28.075*mm;
G4Paraboloid *solidColliA2B2= new G4Paraboloid("solidColliA2B2",ColliA2B2_DemiZ,ColliA2B2_Rayon1,ColliA2B2_Rayon2);

G4double ColliA2B1_Rayon =70.50*mm;
G4double ColliA2B1_DemiZ =28.75*mm;
G4Tubs *solidColliA2B1= new G4Tubs("solidColliA2B1",0.0*mm ,ColliA2B1_Rayon,ColliA2B1_DemiZ ,Angle_Debut, Angle_Fin);

G4ThreeVector  ColliA2_B2_to_B1_decalage_z(0,0,-1.35*mm);
G4RotationMatrix* ColliA2_B2_to_B1_Rotation_z = new G4RotationMatrix;
ColliA2_B2_to_B1_Rotation_z->rotateX(180.*deg);
G4SubtractionSolid* SolidCollimatorA2 = new G4SubtractionSolid("SolidCollimatorA2",solidColliA2B1,solidColliA2B2,ColliA2_B2_to_B1_Rotation_z,ColliA2_B2_to_B1_decalage_z);
//-------------	LOGICAL VOLUME OF COLLIMATOR_A2  ----------------------------
 CollimatorA2 = new G4LogicalVolume(SolidCollimatorA2, 	 //its solid
			  WNICU, 		 //its material
			"CollimatorA2" ,		 //its name
			 0,0,0);


CollimatorA2->SetVisAttributes(WNICU_attribut);

 //-------------PHYSICAL VOLUME OF COLLIMATOR_A1  ----------------------------

G4double CollimatorA1_Centre_z =dec+44.5*mm;
  pCollimatorA1= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, CollimatorA1_Centre_z)),
		 "pCollimatorA1",   //its name  (2nd constructor)
	CollimatorA1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);                       //copy number 
 //-------------PHYSICAL VOLUME OF COLLIMATOR_A2  ----------------------------

G4double CollimatorA2_Centre_z =CollimatorA1_Centre_z+ColliA1_DemiZ+ColliA2B1_DemiZ;
  pCollimatorA2= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, CollimatorA2_Centre_z)),
		 "pCollimatorA2",   //its name  (2nd constructor)
	CollimatorA2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);   

G4Region *regCollimatorA2;
       regCollimatorA2= new G4Region("CollimatorA2Region");
        G4ProductionCuts* cuts = new G4ProductionCuts;
        cuts->SetProductionCut(0.0185*mm,"gamma");
        cuts->SetProductionCut(0.0412*mm,"e-");
        cuts->SetProductionCut(0.0412*mm,"e+");
        regCollimatorA2->SetProductionCuts(cuts);

     CollimatorA2->SetRegion(regCollimatorA2);
        regCollimatorA2->AddRootLogicalVolume(CollimatorA2);

}

// ==================================================================================
//  THIS METHOD SET-UP THE SECONDARY COLLIMATOR
// ==================================================================================
void volumeConstruction::secondaryCollimatorConstructor()
{

//-------------	ATTRIBUT AND MATERIAL OF SECONDARY COLLIMATOR  ----------------------------

 G4VisAttributes *Pb_attribut=this->AttributName("Pb_attribut");
 G4Material *Pb=this->MaterialName("Pb");

//-------------	SOLID VOLUME OF SECONDARY COLLIMATOR  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix

G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;
G4double ColliSecA1_Rayon1 =44.03*mm;
G4double ColliSecA1_Rayon2 =54.88*mm;
G4double ColliSecA1_DemiZ =20.28*mm;
G4Paraboloid *solidColliSecA1= new G4Paraboloid("solidColliSecA1",ColliSecA1_DemiZ,ColliSecA1_Rayon1,ColliSecA1_Rayon2);

G4double ColliSecA2_Rayon =120.5*mm;
G4double ColliSecA2_DemiZ =20.25*mm;
G4Tubs *solidColliSecA2= new G4Tubs("solidColliSecA2",0.0*mm ,ColliSecA2_Rayon,ColliSecA2_DemiZ ,Angle_Debut, Angle_Fin);
G4ThreeVector  ColliSec_A1_to_A2_decalage_z(0,0,0*mm);
G4RotationMatrix* ColliSec_A1_to_A2_Rotation_z = new G4RotationMatrix;
ColliSec_A1_to_A2_Rotation_z->rotateX(180.*deg);
G4SubtractionSolid* SolidColliSec = new G4SubtractionSolid("ColliSec",solidColliSecA2,solidColliSecA1,ColliSec_A1_to_A2_Rotation_z, ColliSec_A1_to_A2_decalage_z);

//-------------	LOGIC VOLUME OF SECONDARY COLLIMATOR  ----------------------------
ColliSec = new G4LogicalVolume(SolidColliSec, 	 //its solid
			  Pb, 		 //its material
			"ColliSec" ,		 //its name
			 0,0,0);

ColliSec->SetVisAttributes(Pb_attribut);

G4double  ColliSec_Centre_z =dec+181.25*mm;
 pColliSec= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ColliSec_Centre_z)),
		 "pColliSec",   //its name  (2nd constructor)
	ColliSec,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4Region *regColliSec;
   regColliSec= new G4Region("ColliSecRegion");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(0.0212*mm,"gamma");
        cuts1->SetProductionCut(0.0652*mm,"e-");
        cuts1->SetProductionCut(0.0652*mm,"e+");
        regColliSec->SetProductionCuts(cuts1);

    ColliSec->SetRegion(regColliSec);
        regColliSec->AddRootLogicalVolume(ColliSec);
        
}
void volumeConstruction::flatteningFilteConstructor()
{
//-------------	ATTRIBUT AND MATERIAL OF FLATTENING FILTER   ----------------------------
 G4VisAttributes *Stainless_steel_attribut=this->AttributName("Stainless_steel_attribut");
G4Material *Stainless_steel=this->MaterialName("Stainless_steel");
//-------------	SOLID VOLUME OF FLATTENING FILTER  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix

G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;
G4double Cone1_Rayon1_max =27.*mm;
G4double Cone1_Rayon1_min =0*mm;
G4double Cone1_Rayon2_max =38*mm;
G4double Cone1_Rayon2_min =0.0*mm;
G4double Cone1_DemiZ =6.35*mm;

G4Cons* Cone1 = new G4Cons("Cone1",Cone1_Rayon1_min,Cone1_Rayon1_max,Cone1_Rayon2_min,Cone1_Rayon2_max,Cone1_DemiZ ,Angle_Debut, Angle_Fin);

//cone2

G4double Cone2_Rayon1_max =4.*mm;
G4double Cone2_Rayon1_min =0*mm;
G4double Cone2_Rayon2_max =27*mm;
G4double Cone2_Rayon2_min =0*mm;
G4double Cone2_DemiZ =15.3*mm;

G4Cons* Cone2 = new G4Cons("Cone2",Cone2_Rayon1_min,Cone2_Rayon1_max,Cone2_Rayon2_min,Cone2_Rayon2_max, Cone2_DemiZ,Angle_Debut, Angle_Fin);


//cone3

G4double Cone3_Rayon1_max =0.00*mm;
G4double Cone3_Rayon1_min =0*mm;
G4double Cone3_Rayon2_max =4*mm;
G4double Cone3_Rayon2_min =0*mm;
G4double Cone3_DemiZ =1.25*mm;

G4Cons* Cone3 = new G4Cons("Cone3",Cone3_Rayon1_min,Cone3_Rayon1_max,Cone3_Rayon2_min,Cone3_Rayon2_max ,Cone3_DemiZ,Angle_Debut, Angle_Fin);

G4ThreeVector  Cone_2_to_1_decalage_z(0,0,-Cone2_DemiZ-Cone1_DemiZ);
G4RotationMatrix* Cone_2_to_1_Rotation_z = new G4RotationMatrix;
Cone_2_to_1_Rotation_z->rotateX(0.*deg);


G4UnionSolid* Cone1U2 = new G4UnionSolid("Cone1U2",Cone1,Cone2,Cone_2_to_1_Rotation_z,Cone_2_to_1_decalage_z);


G4ThreeVector  Cone_3_to_1U2_decalage_z(0,0,-38.2*mm);
G4RotationMatrix* Cone_3_to_1U2_Rotation_z = new G4RotationMatrix;
Cone_3_to_1U2_Rotation_z->rotateX(0.*deg);

G4UnionSolid* Cone1U2U3 = new G4UnionSolid("Cone1U2U3",Cone1U2,Cone3,Cone_3_to_1U2_Rotation_z,Cone_3_to_1U2_decalage_z);

//cylindre1
//cylindre1A

G4double Cylin1A_Rayon =54.0*mm;
G4double Cylin1A_DemiZ =3.25*mm;


G4Tubs *solidCylin1A= new G4Tubs("solidCylin1A",0.0*mm , Cylin1A_Rayon,Cylin1A_DemiZ ,Angle_Debut, Angle_Fin  );

//cylindre1B

G4double Cylin1B_Rayon =43.0*mm;
G4double Cylin1B_DemiZ =3.25*mm;


G4Tubs *solidCylin1B= new G4Tubs("solidCylin1B",0.0*mm , Cylin1B_Rayon,Cylin1B_DemiZ ,Angle_Debut, Angle_Fin  );


G4ThreeVector Cylin1_B_to_A_decalage_z(0,0,0*mm);
G4RotationMatrix* Cylin1_B_to_A_Rotation_z = new G4RotationMatrix;
Cylin1_B_to_A_Rotation_z->rotateX(0.*deg);
G4SubtractionSolid* SolidCylin1 = new G4SubtractionSolid("SolidCylin1",solidCylin1A,solidCylin1B,Cylin1_B_to_A_Rotation_z,Cylin1_B_to_A_decalage_z);

//cylindre2

G4double Cylin2_Rayon =45.0*mm;
G4double Cylin2_DemiZ =0.5*mm;


G4Tubs *solidCylin2= new G4Tubs("solidCylin2",0.0*mm , Cylin2_Rayon,Cylin2_DemiZ ,Angle_Debut, Angle_Fin  );

//FFILTER	

G4ThreeVector Cylin_1_to_2_decalage_z(0,0,-Cylin1B_DemiZ-Cylin2_DemiZ)
;
G4RotationMatrix* Cylin_1_to_2_Rotation_z = new G4RotationMatrix;
Cylin_1_to_2_Rotation_z->rotateX(0.*deg);

G4UnionSolid* SolidCylin = new G4UnionSolid("SolidCylin",solidCylin2,SolidCylin1,Cylin_1_to_2_Rotation_z,Cylin_1_to_2_decalage_z);


G4double demiZ_cone=Cone1_DemiZ+Cone2_DemiZ+Cone3_DemiZ;
G4ThreeVector Cone1U2U3_to_SolidCylin_decalage_z(0,0,-Cylin2_DemiZ-demiZ_cone-1*mm+ (Cone2_DemiZ));

G4RotationMatrix* Cone1U2U3_to_SolidCylin_Rotation_z = new G4RotationMatrix;
Cone1U2U3_to_SolidCylin_Rotation_z->rotateX(0.*deg);

G4UnionSolid* SolidFFILTER = new G4UnionSolid("SolidFFILTER",SolidCylin,Cone1U2U3,Cone1U2U3_to_SolidCylin_Rotation_z,Cone1U2U3_to_SolidCylin_decalage_z);



//-------------	LOGICAL VOLUME OF FLATTENING FILTER  ----------------------------
Cone = new G4LogicalVolume(SolidFFILTER, 	 //its solid
			  Stainless_steel, 		 //its material
			"Cone" ,		 //its name
			 0,0,0);

Cone->SetVisAttributes(Stainless_steel_attribut);
G4double Cone_Centre_z=dec+150.0*mm;

//-------------	PHISICAL VOLUME OF FLATTENING FILTER  ----------------------------
//ColliSec_Centre_z+(2*demiZ_cone)-Colli_Cone_z+3.975*mm;
  pCone= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, Cone_Centre_z)),
		 "pCone",   //its name  (2nd constructor)
	Cone,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4Region *regCone;
   regCone= new G4Region("ConeRegion");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(0.252*mm,"gamma");
        cuts1->SetProductionCut(0.0675*mm,"e-");
        cuts1->SetProductionCut(0.0675*mm,"e+");
        regCone->SetProductionCuts(cuts1);

    Cone->SetRegion(regCone);
        regCone->AddRootLogicalVolume(Cone);
        


}

void volumeConstruction::ionizationChamberConstructor()
{
//-------------	ATTRIBUT AND MATERIAL OF IONIZATION CHAMBER  ----------------------------

 G4VisAttributes *Al_attribut=this->AttributName("Al_attribut");
 G4Material *Al=this->MaterialName("Al");
 G4VisAttributes *Kapton_attribut=this->AttributName("Kapton_attribut");
 G4Material *Kapton=this->MaterialName("Kapton");
//-------------	SOLID VOLUME OF IONIZATION CHAMBER  ----------------------------
G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;
G4double  z_between_ColliSec_and_ChI=1*mm;
G4double  z_between_two_SameWin_ChI=1*mm;
 G4double z_between_two_DiffWin_ChI=2*mm;
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 0.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix
	
//------- 2 X ChI_ElementTypeA:  R = 54 mm aliminum epaisseur 1.6 10-5 cm + kapton epaisseur  25 10-9 m
G4double angestron_to_mm = 0.0000001;
G4double ChI_ElementTypeA_1_Rayon =54*mm;
G4double ChI_ElementTypeA_1_DemiZ= 1600*angestron_to_mm*mm;
G4Tubs *solidChI_ElementTypeA_1= new G4Tubs("solidChI_ElementTypeA_1",0.0*mm , ChI_ElementTypeA_1_Rayon,ChI_ElementTypeA_1_DemiZ ,Angle_Debut, Angle_Fin  );
G4double ChI_ElementTypeA_2_Rayon =54*mm;
G4double ChI_ElementTypeA_2_DemiZ= 250*angestron_to_mm*mm;
G4Tubs *solidChI_ElementTypeA_2= new G4Tubs("solidChI_ElementTypeA_2",0.0*mm ,ChI_ElementTypeA_2_Rayon,ChI_ElementTypeA_2_DemiZ ,Angle_Debut, Angle_Fin  );

G4ThreeVector ChI_ElementTypeA_1_to_2_decalage_z(0,0,-ChI_ElementTypeA_1_Rayon-ChI_ElementTypeA_2_Rayon);
G4RotationMatrix* ChI_ElementTypeA_1_to_2_Rotation_z = new G4RotationMatrix;
ChI_ElementTypeA_1_to_2_Rotation_z->rotateX(0.*deg);

LoGicChI_ElementTypeA1 = new G4LogicalVolume(solidChI_ElementTypeA_1, 	 //its solid
			 Al, 		 //its material
			"LoGicChI_ElementTypeA1" ,		 //its name
			 0,0,0);
LoGicChI_ElementTypeA1->SetVisAttributes(Al_attribut);

LoGicChI_ElementTypeA2 = new G4LogicalVolume(solidChI_ElementTypeA_2, 	 //its solid
			Kapton, 		 //its material
			"LoGicChI_ElementTypeA2" ,		 //its name
		0,0,0);
LoGicChI_ElementTypeA2->SetVisAttributes(Kapton_attribut);

//------- 4 X ChI_ElementTypeB:  R = 54 mm aliminum epaisseur 8 10-6 cm + kapton epaisseur  25 10-9 m

G4double ChI_ElementTypeB_1_Rayon =54*mm;
G4double ChI_ElementTypeB_1_DemiZ = 800*angestron_to_mm*mm;
G4Tubs *solidChI_ElementTypeB_1= new G4Tubs("solidChI_ElementTypeB_1",0.0*mm , ChI_ElementTypeB_1_Rayon,ChI_ElementTypeB_1_DemiZ ,Angle_Debut, Angle_Fin  );
G4double ChI_ElementTypeB_2_Rayon =54*mm;
G4double ChI_ElementTypeB_2_DemiZ =250*angestron_to_mm*mm;
G4Tubs *solidChI_ElementTypeB_2= new G4Tubs("solidChI_ElementTypeB_2",0.0*mm ,ChI_ElementTypeB_2_Rayon,ChI_ElementTypeB_2_DemiZ ,Angle_Debut, Angle_Fin  );

G4ThreeVector ChI_ElementTypeB_1_to_2_decalage_z(0,0,-ChI_ElementTypeB_1_Rayon-ChI_ElementTypeB_2_Rayon);
G4RotationMatrix* ChI_ElementTypeB_1_to_2_Rotation_z = new G4RotationMatrix;
ChI_ElementTypeB_1_to_2_Rotation_z->rotateX(0.*deg);


 LoGicChI_ElementTypeB1 = new G4LogicalVolume(solidChI_ElementTypeB_1, 	 //its solid
			  Al, 		 //its material
			"LoGicChI_ElementTypeB1" ,		 //its name
			 0,0,0);
LoGicChI_ElementTypeB1->SetVisAttributes(Al_attribut);

LoGicChI_ElementTypeB2 = new G4LogicalVolume(solidChI_ElementTypeB_2, 	 //its solid
			  Kapton, 		 //its material
			"LoGicChI_ElementTypeB2" ,		 //its name
			 0,0,0);
LoGicChI_ElementTypeB2->SetVisAttributes(Kapton_attribut);


//-------------	PHYSICAL VOLUME OF IONIZATION CHAMBER  ----------------------------
G4double ColliSecA1_DemiZ =20.28*mm;

//CHAMBRE D IONISATION
// LAYER1
  G4double ColliSec_Centre_z =dec+181.25*mm;

G4double ChI_ElementTypeB1_Centre_z_layer1=ColliSec_Centre_z+ColliSecA1_DemiZ+z_between_ColliSec_and_ChI
+ChI_ElementTypeB_1_DemiZ;

 pChI_ElementTypeB1_layer1= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB1_Centre_z_layer1)),
		 "pLoGicChI_ElementTypeB1_layer1",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double ChI_ElementTypeB2_Centre_z_layer1=ChI_ElementTypeB1_Centre_z_layer1+ChI_ElementTypeB_2_DemiZ
+ChI_ElementTypeB_1_DemiZ;

 pChI_ElementTypeB2_layer1= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB2_Centre_z_layer1)),
		 "pLoGicChI_ElementTypeB2_layer1",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

// LAYER2

G4double ChI_ElementTypeA1_Centre_z_layer2=ChI_ElementTypeB2_Centre_z_layer1+z_between_two_DiffWin_ChI
+ChI_ElementTypeB_2_DemiZ-ChI_ElementTypeA_1_DemiZ;//

  pChI_ElementTypeA1_layer2= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeA1_Centre_z_layer2)),
		 "pLoGicChI_ElementTypeA1_layer2",   //its name  (2nd constructor)
	LoGicChI_ElementTypeA1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double ChI_ElementTypeA2_Centre_z_layer2=ChI_ElementTypeA1_Centre_z_layer2+ChI_ElementTypeA_2_DemiZ
+ChI_ElementTypeA_1_DemiZ;

  pChI_ElementTypeA2_layer2= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeA2_Centre_z_layer2)),
		 "pLoGicChI_ElementTypeA2_layer2",   //its name  (2nd constructor)
	LoGicChI_ElementTypeA2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

// LAYER3

G4double ChI_ElementTypeB1_Centre_z_layer3=ChI_ElementTypeA_2_DemiZ+ChI_ElementTypeA2_Centre_z_layer2
+z_between_two_DiffWin_ChI+ChI_ElementTypeB_1_DemiZ;//

 pChI_ElementTypeB1_layer3= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB1_Centre_z_layer3)),
		 "pLoGicChI_ElementTypeB1_layer3",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double ChI_ElementTypeB2_Centre_z_layer3=ChI_ElementTypeB1_Centre_z_layer3+ChI_ElementTypeB_2_DemiZ
+ChI_ElementTypeB_1_DemiZ;

 pChI_ElementTypeB2_layer3= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB2_Centre_z_layer3)),
		 "pLoGicChI_ElementTypeB2_layer3",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

// LAYER4

G4double ChI_ElementTypeB1_Centre_z_layer4=ChI_ElementTypeB2_Centre_z_layer3+z_between_two_SameWin_ChI
+ChI_ElementTypeB_2_DemiZ-ChI_ElementTypeB_1_DemiZ;//

 pChI_ElementTypeB1_layer4= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB1_Centre_z_layer4)),
		 "pLoGicChI_ElementTypeB1_layer4",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4Region *regChI_ElementTypeB1;
   regChI_ElementTypeB1= new G4Region("ChI_ElementTypeB1Region");
        G4ProductionCuts* cuts3 = new G4ProductionCuts;
        cuts3->SetProductionCut(2.2*mm,"gamma");
        cuts3->SetProductionCut(0.172*mm,"e-");
        cuts3->SetProductionCut(0.172*mm,"e+");
        regChI_ElementTypeB1->SetProductionCuts(cuts3);

    LoGicChI_ElementTypeB1->SetRegion(regChI_ElementTypeB1);
        regChI_ElementTypeB1->AddRootLogicalVolume(LoGicChI_ElementTypeB1);

G4double ChI_ElementTypeB2_Centre_z_layer4=ChI_ElementTypeB1_Centre_z_layer4+ChI_ElementTypeB_2_DemiZ
+ChI_ElementTypeB_1_DemiZ;

 pChI_ElementTypeB2_layer4= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB2_Centre_z_layer4)),
		 "pLoGicChI_ElementTypeB2_layer4",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);


// LAYER5

G4double ChI_ElementTypeA1_Centre_z_layer5=ChI_ElementTypeB2_Centre_z_layer4+z_between_two_DiffWin_ChI
+ChI_ElementTypeB_2_DemiZ-ChI_ElementTypeA_1_DemiZ;//

  pChI_ElementTypeA1_layer5= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeA1_Centre_z_layer5)),
		 "pLoGicChI_ElementTypeA1_layer5",   //its name  (2nd constructor)
	LoGicChI_ElementTypeA1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4Region *regChI_ElementTypeA1;
   regChI_ElementTypeA1= new G4Region("ChI_ElementTypeA1Region");
        G4ProductionCuts* cuts5 = new G4ProductionCuts;
        cuts5->SetProductionCut(2.2*mm,"gamma");
        cuts5->SetProductionCut(0.172*mm,"e-");
        cuts5->SetProductionCut(0.172*mm,"e+");
        regChI_ElementTypeA1->SetProductionCuts(cuts5);

    LoGicChI_ElementTypeA1->SetRegion(regChI_ElementTypeA1);
        regChI_ElementTypeA1->AddRootLogicalVolume(LoGicChI_ElementTypeA1);
        
G4double ChI_ElementTypeA2_Centre_z_layer5=ChI_ElementTypeA1_Centre_z_layer5+ChI_ElementTypeA_2_DemiZ
+ChI_ElementTypeA_1_DemiZ;

 pChI_ElementTypeA2_layer5= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeA2_Centre_z_layer5)),
		 "pLoGicChI_ElementTypeA2_layer5",   //its name  (2nd constructor)
	LoGicChI_ElementTypeA2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4Region *regChI_ElementTypeA2;
   regChI_ElementTypeA2= new G4Region("ChI_ElementTypeA2Region");
        G4ProductionCuts* cuts = new G4ProductionCuts;
        cuts->SetProductionCut(21.6*mm,"gamma");
        cuts->SetProductionCut(0.28*mm,"e-");
        cuts->SetProductionCut(0.28*mm,"e+");
        regChI_ElementTypeA2->SetProductionCuts(cuts);

    LoGicChI_ElementTypeA2->SetRegion(regChI_ElementTypeA2);
        regChI_ElementTypeA2->AddRootLogicalVolume(LoGicChI_ElementTypeA2);
//LAYER6

G4double ChI_ElementTypeB1_Centre_z_layer6=ChI_ElementTypeA_2_DemiZ+ChI_ElementTypeA2_Centre_z_layer5
+z_between_two_DiffWin_ChI+ChI_ElementTypeB_1_DemiZ;//

 pChI_ElementTypeB1_layer6= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB1_Centre_z_layer6)),
		 "pLoGicChI_ElementTypeB1_layer6",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double ChI_ElementTypeB2_Centre_z_layer6=ChI_ElementTypeB1_Centre_z_layer6+ChI_ElementTypeB_2_DemiZ
+ChI_ElementTypeB_1_DemiZ;

pChI_ElementTypeB2_layer6= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, ChI_ElementTypeB2_Centre_z_layer6)),
		 "pLoGicChI_ElementTypeB2_layer6",   //its name  (2nd constructor)
	LoGicChI_ElementTypeB2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
//-------------------------------------------------------

G4Region *regChI_ElementTypeB2;
   regChI_ElementTypeB2= new G4Region("ChI_ElementTypeB2Region");
        G4ProductionCuts* cuts1 = new G4ProductionCuts;
        cuts1->SetProductionCut(21.6*mm,"gamma");
        cuts1->SetProductionCut(0.28*mm,"e-");
        cuts1->SetProductionCut(0.28*mm,"e+");
        regChI_ElementTypeB2->SetProductionCuts(cuts1);

    LoGicChI_ElementTypeB2->SetRegion(regChI_ElementTypeB2);
        regChI_ElementTypeB2->AddRootLogicalVolume(LoGicChI_ElementTypeB2);
        


}
void volumeConstruction::AlPlaqueConstructor()
{

//-------------	ATTRIBUT AND MATERIAL OF ALIMUNIUM PLAQUE  ----------------------------

 G4VisAttributes *Al_attribut=this->AttributName("Al_attribut");
 G4Material *Al=this->MaterialName("Al");
//-------------	SOLID VOLUME OF ALIMUNIUM PLAQUE  ----------------------------
G4RotationMatrix rotMatrixpTarget;   // unit rotation matrix
G4double anglepTarget= 180.0*deg;   // rotational angle
rotMatrixpTarget.rotateX(anglepTarget);  // rot matrix
G4double Angle_Debut =0.0*rad;
G4double Angle_Fin =360.0*rad;
G4double PlqAl_Rayon =58.265*mm;
 G4double PlqAl_DemiZ =1.0*mm;


G4Tubs *solidPlqAl= new G4Tubs("solidPlqAl",0.0*mm , PlqAl_Rayon,PlqAl_DemiZ ,Angle_Debut, Angle_Fin  );

//-------------	LOGICAL VOLUME OF ALIMUNIUM PLAQUE  ----------------------------
 LoGiCPlqAl = new G4LogicalVolume(solidPlqAl, 	 //its solid
			  Al, 		 //its material
			"LoGiCPlqAl" ,		 //its name
			 0,0,0);

LoGiCPlqAl->SetVisAttributes(Al_attribut);
 

//-------------	PHYSICAL VOLUME OF ALIMUNIUM PLAQUE  ----------------------------
G4double z_Centre_ChI_PlAl=dec+216.75*mm;


 pLoGiCPlqAl= new G4PVPlacement(G4Transform3D(rotMatrixpTarget,	//rotation 
		 G4ThreeVector(0*cm, 0.0*m, z_Centre_ChI_PlAl)),
		 "pLoGiCPlqAl",   //its name  (2nd constructor)
	LoGiCPlqAl,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4Region *regLoGiCPlqAl;
   regLoGiCPlqAl= new G4Region("LoGiCPlqAlRegion");
        G4ProductionCuts* cuts4 = new G4ProductionCuts;
        cuts4->SetProductionCut(2.2*mm,"gamma");
        cuts4->SetProductionCut(0.172*mm,"e-");
        cuts4->SetProductionCut(0.172*mm,"e+");
        regLoGiCPlqAl->SetProductionCuts(cuts4);

    LoGiCPlqAl->SetRegion(regLoGiCPlqAl);
        regLoGiCPlqAl->AddRootLogicalVolume(LoGiCPlqAl);


}
void volumeConstruction::jawsConstructor()
{
//-------------	ATTRIBUT AND MATERIAL OF JAWS  ----------------------------

 G4VisAttributes *XC10_attribut=this->AttributName("XC10_attribut");
 G4Material *XC10=this->MaterialName("XC10");
 G4VisAttributes *WNICU_attribut=this->AttributName("WNICU_attribut");
 G4Material *WNICU=this->MaterialName("WNICU");
 G4VisAttributes *Pb_attribut=this->AttributName("Pb_attribut");
 G4Material *Pb=this->MaterialName("Pb");


 G4double field=5.*cm;
G4double theta = std::atan(field/(100.*cm));// in rad
//G4double theta = 0;// in rad

 //-------------	SOLID & LOGICAL VOLUME OF JAWS  ----------------------------
//jaws
//jaws X1 & X2

//JAWS COLIMATORS
//JAWS X
//JAW X1 and X2 :6 sub collimator

//----------------Jaw X1_1---------------

  G4double JawX1_1_Dim_y = 101*mm;
  G4double JawX1_1_Dim_x = 50.5*mm;
  G4double JawX1_1_Dim_z = 1.5*mm;
  G4Box* JawX1_1 = new G4Box("JawX1_1",JawX1_1_Dim_x,JawX1_1_Dim_y,JawX1_1_Dim_z);
  LoGiCJawX1_1 = new G4LogicalVolume(JawX1_1,XC10,"LoGiCJawX1_1",0,0,0);
  LoGiCJawX1_1->SetVisAttributes(XC10_attribut);
G4Region *regLoGiCJawX1_1;
   regLoGiCJawX1_1= new G4Region("LoGiCJawX1_1Region");
        G4ProductionCuts* cuts31= new G4ProductionCuts;
        cuts31->SetProductionCut(0.25*mm,"gamma");
        cuts31->SetProductionCut(0.068*mm,"e-");
        cuts31->SetProductionCut(0.068*mm,"e+");
        regLoGiCJawX1_1->SetProductionCuts(cuts31);
        LoGiCJawX1_1->SetRegion(regLoGiCJawX1_1);
        regLoGiCJawX1_1->AddRootLogicalVolume(LoGiCJawX1_1);
//----------------Jaw X1_2---------------

  G4double JawX1_2_Dim_y = 101*mm;
  G4double JawX1_2_Dim_x = 50.5*mm;
  G4double JawX1_2_Dim_z = 17.5*mm;
  G4Box* JawX1_2 = new G4Box("JawX1_2",JawX1_2_Dim_x,JawX1_2_Dim_y,JawX1_2_Dim_z);
  LoGiCJawX1_2 = new G4LogicalVolume(JawX1_2,Pb,"LoGiCJawX1_2",0,0,0);
  LoGiCJawX1_2->SetVisAttributes(Pb_attribut);
G4Region *regLoGiCJawX1_2;
   regLoGiCJawX1_2= new G4Region("LoGiCJawX1_2Region");
        G4ProductionCuts* cuts30= new G4ProductionCuts;
        cuts30->SetProductionCut(0.0212*mm,"gamma");
        cuts30->SetProductionCut(0.0652*mm,"e-");
        cuts30->SetProductionCut(0.0652*mm,"e+");
           regLoGiCJawX1_2->SetProductionCuts(cuts30);
     LoGiCJawX1_2->SetRegion(regLoGiCJawX1_2);
           regLoGiCJawX1_2->AddRootLogicalVolume(LoGiCJawX1_2);

//----------------Jaw X1_3---------------

  G4double JawX1_3_Dim_y = 101*mm;
  G4double JawX1_3_Dim_x = 50.5*mm;
  G4double JawX1_3_Dim_z = 17.5*mm;
  G4Box* JawX1_3 = new G4Box("JawX1_3",JawX1_3_Dim_x,JawX1_3_Dim_y,JawX1_3_Dim_z);
  LoGiCJawX1_3 = new G4LogicalVolume(JawX1_3, WNICU,"LoGiCJawX1_3",0,0,0);
  LoGiCJawX1_3->SetVisAttributes(WNICU_attribut);
G4Region *regLoGiCJawX1_3;
   regLoGiCJawX1_3= new G4Region("LoGiCJawX1_3Region");
        G4ProductionCuts* cuts10= new G4ProductionCuts;
        cuts10->SetProductionCut(0.0185*mm,"gamma");
        cuts10->SetProductionCut(0.0412*mm,"e-");
        cuts10->SetProductionCut(0.0412*mm,"e+");
           regLoGiCJawX1_3->SetProductionCuts(cuts10);
     LoGiCJawX1_3->SetRegion(regLoGiCJawX1_3);
           regLoGiCJawX1_3->AddRootLogicalVolume(LoGiCJawX1_3);
//----------------Jaw X1_4---------------

  G4double JawX1_4_Dim_y = 101*mm;
  G4double JawX1_4_Dim_x = 50.5*mm;
  G4double JawX1_4_Dim_z = 13.5*mm;
  G4Box* JawX1_4 = new G4Box("JawX1_4",JawX1_4_Dim_x,JawX1_4_Dim_y,JawX1_4_Dim_z);
  LoGiCJawX1_4 = new G4LogicalVolume(JawX1_4,Pb,"LoGiCJawX1_4",0,0,0);
  LoGiCJawX1_4->SetVisAttributes(Pb_attribut);
G4Region *regLoGiCJawX1_4;
   regLoGiCJawX1_4= new G4Region("LoGiCJawX1_4Region");
        G4ProductionCuts* cuts17154= new G4ProductionCuts;
        cuts17154->SetProductionCut(0.0212*mm,"gamma");
        cuts17154->SetProductionCut(0.0652*mm,"e-");
        cuts17154->SetProductionCut(0.0652*mm,"e+");
           regLoGiCJawX1_4->SetProductionCuts(cuts17154);
     LoGiCJawX1_4->SetRegion(regLoGiCJawX1_4);
           regLoGiCJawX1_4->AddRootLogicalVolume(LoGiCJawX1_4);
//----------------Jaw X1_5---------------

  G4double JawX1_5_Dim_y = 101*mm;
  G4double JawX1_5_Dim_x = 50.5*mm;
  G4double JawX1_5_Dim_z = 5*mm;
  G4Box* JawX1_5 = new G4Box("JawX1_5",JawX1_5_Dim_x, JawX1_5_Dim_y,JawX1_5_Dim_z);
  LoGiCJawX1_5 = new G4LogicalVolume(JawX1_5, Pb,"LoGiCJawX1_5",0,0,0);
  LoGiCJawX1_5->SetVisAttributes(Pb_attribut);
G4Region *regLoGiCJawX1_5;
   regLoGiCJawX1_5= new G4Region("LoGiCJawX1_5Region");
        G4ProductionCuts* cuts1715= new G4ProductionCuts;
        cuts1715->SetProductionCut(0.0212*mm,"gamma");
        cuts1715->SetProductionCut(0.0652*mm,"e-");
        cuts1715->SetProductionCut(0.0652*mm,"e+");
           regLoGiCJawX1_5->SetProductionCuts(cuts1715);
     LoGiCJawX1_5->SetRegion(regLoGiCJawX1_5);
           regLoGiCJawX1_5->AddRootLogicalVolume(LoGiCJawX1_5);
//----------------Jaw X1_6---------------

  G4double JawX1_6_Dim_y = 101*mm;
  G4double JawX1_6_Dim_x = 50.5*mm;
  G4double JawX1_6_Dim_z = 2.5*mm;
  G4Box* JawX1_6 = new G4Box("JawX1_6",JawX1_6_Dim_x,JawX1_6_Dim_y,JawX1_6_Dim_z);
  LoGiCJawX1_6 = new G4LogicalVolume(JawX1_6,WNICU,"LoGiCJawX1_6",0,0,0);
  LoGiCJawX1_6->SetVisAttributes(WNICU_attribut);
G4Region *regLoGiCJawX1_6;
   regLoGiCJawX1_6= new G4Region("LoGiCJawX1_6Region");
        G4ProductionCuts* cuts171= new G4ProductionCuts;
        cuts171->SetProductionCut(0.0185*mm,"gamma");
        cuts171->SetProductionCut(0.0412*mm,"e-");
        cuts171->SetProductionCut(0.0412*mm,"e+");
           regLoGiCJawX1_6->SetProductionCuts(cuts171);
     LoGiCJawX1_6->SetRegion(regLoGiCJawX1_6);
           regLoGiCJawX1_6->AddRootLogicalVolume(LoGiCJawX1_6);

//JAWS Y
//JAW Y1 :4 sub collimator

//----------------Jaw Y1_1---------------

  G4double JawY1_1_Dim_x = 101*mm;
  G4double JawY1_1_Dim_y = 50.5*mm;
  G4double JawY1_1_Dim_z = 15.5*mm;
  G4Box* JawY1_1 = new G4Box("JawX1_1",JawY1_1_Dim_x, JawY1_1_Dim_y,JawY1_1_Dim_z);
  LoGiCJawY1_1 = new G4LogicalVolume(JawY1_1,Pb,"LoGiCJawY1_1",0,0,0);
  LoGiCJawY1_1->SetVisAttributes(Pb_attribut);
G4Region *regLoGiCJawY1_1;
   regLoGiCJawY1_1= new G4Region("LoGiCJawY1_1Region");
        G4ProductionCuts* cuts6 = new G4ProductionCuts;
        cuts6->SetProductionCut(0.0212*mm,"gamma");
        cuts6->SetProductionCut(0.0652*mm,"e-");
        cuts6->SetProductionCut(0.0652*mm,"e+");
           regLoGiCJawY1_1->SetProductionCuts(cuts6);

      LoGiCJawY1_1->SetRegion(   regLoGiCJawY1_1);
           regLoGiCJawY1_1->AddRootLogicalVolume(LoGiCJawY1_1);

//----------------Jaw Y1_2---------------

  G4double JawY1_2_Dim_x = 101*mm;
  G4double JawY1_2_Dim_y = 50.5*mm;
  G4double JawY1_2_Dim_z = 17.5*mm;
  G4Box* JawY1_2 = new G4Box("JawY1_2",JawY1_2_Dim_x,JawY1_2_Dim_y,JawY1_2_Dim_z);
  LoGiCJawY1_2 = new G4LogicalVolume(JawY1_2,WNICU,"LoGiCJawY1_2",0,0,0);
  LoGiCJawY1_2->SetVisAttributes(WNICU_attribut);
G4Region *regLoGiCJawY1_2;
   regLoGiCJawY1_2= new G4Region("LoGiCJawY1_2Region");
        G4ProductionCuts* cuts7 = new G4ProductionCuts;
        cuts7->SetProductionCut(0.0185*mm,"gamma");
        cuts7->SetProductionCut(0.0412*mm,"e-");
        cuts7->SetProductionCut(0.0412*mm,"e+");
           regLoGiCJawY1_2->SetProductionCuts(cuts7);

      LoGiCJawY1_2->SetRegion(   regLoGiCJawY1_2);
           regLoGiCJawY1_2->AddRootLogicalVolume(LoGiCJawY1_2);

//----------------Jaw Y1_3---------------

  G4double JawY1_3_Dim_x = 101*mm;
  G4double JawY1_3_Dim_y = 50.5*mm;
  G4double JawY1_3_Dim_z = 10.5*mm;
  G4Box* JawY1_3 = new G4Box("JawY1_3",JawY1_3_Dim_x,JawY1_3_Dim_y,JawY1_3_Dim_z);
  LoGiCJawY1_3 = new G4LogicalVolume(JawY1_3,Pb,"LoGiCJawY1_3",0,0,0);
  LoGiCJawY1_3->SetVisAttributes(Pb_attribut);
G4Region *regLoGiCJawY1_3;
   regLoGiCJawY1_3= new G4Region("LoGiCJawY1_3Region");
        G4ProductionCuts* cuts37 = new G4ProductionCuts;
        cuts37->SetProductionCut(0.0212*mm,"gamma");
        cuts37->SetProductionCut(0.0652*mm,"e-");
        cuts37->SetProductionCut(0.0652*mm,"e+");
           regLoGiCJawY1_3->SetProductionCuts(cuts37);

      LoGiCJawY1_3->SetRegion(   regLoGiCJawY1_3);
           regLoGiCJawY1_3->AddRootLogicalVolume(LoGiCJawY1_3);
//----------------Jaw Y1_4---------------

  G4double JawY1_4_Dim_x = 101*mm;
  G4double JawY1_4_Dim_y = 50.5*mm;
  G4double JawY1_4_Dim_z = 7.5*mm;
  G4Box* JawY1_4 = new G4Box("JawY1_4",JawY1_4_Dim_x,JawY1_4_Dim_y,JawY1_4_Dim_z);
  LoGiCJawY1_4 = new G4LogicalVolume(JawY1_4,WNICU,"LoGiCJawY1_4",0,0,0);
  LoGiCJawY1_4->SetVisAttributes( WNICU_attribut);
G4Region *regLoGiCJawY1_4;
   regLoGiCJawY1_4= new G4Region("LoGiCJawY1_4Region");
        G4ProductionCuts* cuts71= new G4ProductionCuts;
        cuts71->SetProductionCut(0.0185*mm,"gamma");
        cuts71->SetProductionCut(0.0412*mm,"e-");
        cuts71->SetProductionCut(0.0412*mm,"e+");
           regLoGiCJawY1_4->SetProductionCuts(cuts71);

      LoGiCJawY1_4->SetRegion(   regLoGiCJawY1_4);
           regLoGiCJawY1_4->AddRootLogicalVolume(LoGiCJawY1_4);

 //-------------	PHYSICAL VOLUME OF JAWS  ----------------------------

G4double JawX1_1_centre_y=0*mm;
G4double distance_Source_JawX1_1=275.5*mm+JawX1_1_Dim_z;
G4double dx_x1_1=distance_Source_JawX1_1*tan(theta);
G4double JawX1_1_centre_x=JawX1_1_Dim_x+dx_x1_1;
G4double JawX1_1_centre_z=275.5*mm+JawX1_1_Dim_z+dec;
G4RotationMatrix RMX1_1,RMX2_1;   // unit rotation matrix
G4double angleX1_1 = -theta;   // rotational angle
RMX1_1.rotateY(-angleX1_1);  // rot matrix
RMX2_1.rotateY(angleX1_1);  // rot matrix
 pLoGiCJawX1_1= new G4PVPlacement(G4Transform3D(RMX1_1,	//rotation 
		 G4ThreeVector(JawX1_1_centre_x, JawX1_1_centre_y, JawX1_1_centre_z)),
		 "pLoGiCJawX1_1",   //its name  (2nd constructor)
	LoGiCJawX1_1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
 pLoGiCJawX2_1= new G4PVPlacement(G4Transform3D(RMX2_1,	//rotation 
		 G4ThreeVector(-JawX1_1_centre_x, JawX1_1_centre_y, JawX1_1_centre_z)),
		 "pLoGiCJawX2_1",   //its name  (2nd constructor)
	LoGiCJawX1_1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4double distance_Source_JawX1_2=278.5*mm+JawX1_2_Dim_z;
G4double dx_x1_2=distance_Source_JawX1_2*tan(theta);
G4RotationMatrix RMX1_2,RMX2_2;   // unit rotation matrix
G4double angleX1_2 = theta;   // rotational angle
RMX1_2.rotateY(angleX1_2);  // rot matrix
RMX2_2.rotateY(-angleX1_2);  // rot matrix
G4double JawX1_2_centre_x=dx_x1_2+JawX1_2_Dim_x;
G4double JawX1_2_centre_y=0*mm;
G4double JawX1_2_centre_z= 278.5*mm+JawX1_2_Dim_z+dec;

 pLoGiCJawX1_2= new G4PVPlacement(G4Transform3D(RMX1_2,	//rotation 
		 G4ThreeVector(JawX1_2_centre_x, JawX1_2_centre_y, JawX1_2_centre_z)),
		 "pLoGiCJawX1_2",   //its name  (2nd constructor)
	LoGiCJawX1_2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

pLoGiCJawX2_2= new G4PVPlacement(G4Transform3D(RMX2_2,	//rotation 
		 G4ThreeVector(-JawX1_2_centre_x, JawX1_2_centre_y, JawX1_2_centre_z)),
		 "pLoGiCJawX2_2",   //its name  (2nd constructor)
	LoGiCJawX1_2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double distance_Source_JawX1_3=313.5*mm+JawX1_3_Dim_z;
G4double dx_x1_3=distance_Source_JawX1_3*tan(theta);
G4double JawX1_3_centre_x=dx_x1_3+JawX1_3_Dim_x;
G4double JawX1_3_centre_y=0*mm;

G4double JawX1_3_centre_z=313.5*mm+JawX1_3_Dim_z+dec;
G4RotationMatrix RMX1_3,RMX2_3;   // unit rotation matrix
G4double angleX1_3 = theta;   // rotational angle
RMX1_3.rotateY(angleX1_3);
RMX2_3.rotateY(-angleX1_3);
 pLoGiCJawX1_3= new G4PVPlacement(G4Transform3D(RMX1_3,	//rotation 
		 G4ThreeVector(JawX1_3_centre_x, JawX1_3_centre_y, JawX1_3_centre_z)),
		 "pLoGiCJawX1_3",   //its name  (2nd constructor)
	LoGiCJawX1_3,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

pLoGiCJawX2_3= new G4PVPlacement(G4Transform3D(RMX2_3,	//rotation 
		 G4ThreeVector(-JawX1_3_centre_x, JawX1_3_centre_y, JawX1_3_centre_z)),
		 "pLoGiCJawX2_3",   //its name  (2nd constructor)
	LoGiCJawX1_3,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);




G4double distance_Source_JawX1_4=348.5*mm+JawX1_4_Dim_z;
G4double dx_x1_4=distance_Source_JawX1_4*tan(theta);

G4RotationMatrix RMX1_4,RMX2_4;   // unit rotation matrix
G4double angleX1_4 = theta;   // rotational angle
RMX1_4.rotateY(angleX1_4);
RMX2_4.rotateY(-angleX1_4);
G4double JawX1_4_centre_x=JawX1_4_Dim_x+dx_x1_4;
G4double JawX1_4_centre_y=0*mm;
G4double JawX1_4_centre_z=348.5*mm+JawX1_4_Dim_z+dec;

 pLoGiCJawX1_4= new G4PVPlacement(G4Transform3D(RMX1_4,	//rotation 
		 G4ThreeVector(JawX1_4_centre_x, JawX1_4_centre_y, JawX1_4_centre_z)),
		 "pLoGiCJawX1_4",   //its name  (2nd constructor)
	LoGiCJawX1_4,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

pLoGiCJawX2_4= new G4PVPlacement(G4Transform3D(RMX2_4,	//rotation 
		 G4ThreeVector(-JawX1_4_centre_x, JawX1_4_centre_y, JawX1_4_centre_z)),
		 "pLoGiCJawX2_4",   //its name  (2nd constructor)
	LoGiCJawX1_4,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
G4double distance_Source_JawX1_5=472.5*mm+JawX1_5_Dim_z;
G4double dx_x1_5=distance_Source_JawX1_5*tan(theta);
G4RotationMatrix RMX1_5,RMX2_5;   // unit rotation matrix
G4double angleX1_5 = theta;   // rotational angle
RMX1_5.rotateY(angleX1_5);
RMX2_5.rotateY(-angleX1_5);
G4double JawX1_5_centre_x=JawX1_5_Dim_x+dx_x1_5;
G4double JawX1_5_centre_y=0*mm;
G4double JawX1_5_centre_z=472.5*mm+JawX1_5_Dim_z+dec;

 pLoGiCJawX1_5= new G4PVPlacement(G4Transform3D(RMX1_5,	//rotation 
		 G4ThreeVector(JawX1_5_centre_x, JawX1_5_centre_y, JawX1_5_centre_z)),   //its name  (2nd c
	 "pLoGiCJawX1_5",LoGiCJawX1_5,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

pLoGiCJawX2_5= new G4PVPlacement(G4Transform3D(RMX2_5,	//rotation 
		 G4ThreeVector(-JawX1_5_centre_x, JawX1_5_centre_y, JawX1_5_centre_z)),   //its name  (2nd c
	 "pLoGiCJawX2_5",LoGiCJawX1_5,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

G4double distance_Source_JawX1_6=482.5*mm+JawX1_6_Dim_z;

  G4double dx_x1_6=distance_Source_JawX1_6*tan(theta);
  G4RotationMatrix RMX1_6,RMX2_6;   // unit rotation matrix
  G4double angleX1_6 = theta;   // rotational angle
  RMX1_6.rotateY(angleX1_6);
  RMX2_6.rotateY(-angleX1_6);
  G4double JawX1_6_centre_x=JawX1_6_Dim_x+dx_x1_6;
  G4double JawX1_6_centre_y=0*mm;
G4double JawX1_6_centre_z=482.5*mm+JawX1_6_Dim_z+dec;

   pLoGiCJawX1_6= new G4PVPlacement(G4Transform3D(RMX1_6,	//rotation 
		 G4ThreeVector(JawX1_6_centre_x, JawX1_6_centre_y, JawX1_6_centre_z)),
		 "pLoGiCJawX1_6",   //its name  (2nd constructor)
	LoGiCJawX1_6,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

   pLoGiCJawX2_6= new G4PVPlacement(G4Transform3D(RMX2_6,	//rotation 
		 G4ThreeVector(-JawX1_6_centre_x, JawX1_6_centre_y, JawX1_6_centre_z)),
		 "pLoGiCJawX2_6",   //its name  (2nd constructor)
	LoGiCJawX1_6,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
//jaws Y
   G4double distance_Source_JawY1_1=380.5*mm+JawY1_1_Dim_z;
   G4double dy_y1_1=distance_Source_JawY1_1*tan(theta);
   G4RotationMatrix RMY1_1,RMY2_1;   // unit rotation matrix
   G4double angleY1_1 = -theta;   // rotational angle
   RMY1_1.rotateX(angleY1_1);
   RMY2_1.rotateX(-angleY1_1);
   G4double JawY1_1_centre_x=0*mm;
   G4double jawY1_1_centre_y=JawY1_1_Dim_y+dy_y1_1;
   //G4double jawY1_1_centre_z=JawX1_4_centre_z-JawY1_1_Dim_z-distance_JawX1_4_JawY1_1-JawX1_4_Dim_z;
G4double jawY1_1_centre_z=380.5*mm+JawY1_1_Dim_z+dec;
   pLoGiCJawY1_1= new G4PVPlacement(G4Transform3D( RMY1_1,	//rotation 
		 G4ThreeVector(JawY1_1_centre_x, jawY1_1_centre_y, jawY1_1_centre_z)),
		 "pLoGiCJawY1_1",   //its name  (2nd constructor)
	LoGiCJawY1_1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

   pLoGiCJawY2_1= new G4PVPlacement(G4Transform3D( RMY2_1,	//rotation 
		 G4ThreeVector(JawY1_1_centre_x, -jawY1_1_centre_y, jawY1_1_centre_z)),
		 "pLoGiCJawY2_1",   //its name  (2nd constructor)
	LoGiCJawY1_1,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);

   G4double distance_Source_JawY1_2=411.5*mm+JawY1_2_Dim_z;
   G4double dy_y1_2=distance_Source_JawY1_2*tan(theta);
   G4RotationMatrix RMY1_2,RMY2_2;   // unit rotation matrix
   G4double angleY1_2 = -theta;   // rotational angle
   RMY1_2.rotateX(angleY1_2);
   RMY2_2.rotateX(-angleY1_2);
   G4double jawY1_2_centre_x= 0*mm;
   G4double jawY1_2_centre_y=JawY1_2_Dim_y+dy_y1_2;
   G4double jawY1_2_centre_z= 411.5*mm+JawY1_2_Dim_z+dec;

   pLoGiCJawY1_2= new G4PVPlacement(G4Transform3D(RMY1_2,	//rotation 
		 G4ThreeVector(jawY1_2_centre_x, jawY1_2_centre_y, jawY1_2_centre_z)),
		 "pLoGiCJawY1_2",   //its name  (2nd constructor)
	LoGiCJawY1_2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
   pLoGiCJawY2_2= new G4PVPlacement(G4Transform3D(RMY2_2,	//rotation 
		 G4ThreeVector(jawY1_2_centre_x, -jawY1_2_centre_y, jawY1_2_centre_z)),
		 "pLoGiCJawY2_2",   //its name  (2nd constructor)
	LoGiCJawY1_2,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
  G4double distance_Source_JawY1_3=446.5*mm+JawY1_3_Dim_z;
  G4double dy_y1_3=distance_Source_JawY1_3*tan(theta);
  G4RotationMatrix RMY1_3,RMY2_3;   // unit rotation matrix
  G4double angleY1_3 = -theta;   // rotational angle
  RMY1_3.rotateX(angleY1_3);
  RMY2_3.rotateX(-angleY1_3);
  G4double jawY1_3_centre_y=JawY1_3_Dim_y+dy_y1_3;
  G4double jawY1_3_centre_x=0*mm;
  G4double jawY1_3_centre_z=446.5*mm+JawY1_3_Dim_z+dec;

  pLoGiCJawY1_3= new G4PVPlacement(G4Transform3D(RMY1_3,	//rotation 
		 G4ThreeVector(jawY1_3_centre_x, jawY1_3_centre_y, jawY1_3_centre_z)),
		 "pLoGiCJawY1_3",   //its name  (2nd constructor)
	LoGiCJawY1_3,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
  pLoGiCJawY2_3= new G4PVPlacement(G4Transform3D(RMY2_3,	//rotation 
		 G4ThreeVector(jawY1_3_centre_x, -jawY1_3_centre_y, jawY1_3_centre_z)),
		 "pLoGiCJawY2_3",   //its name  (2nd constructor)
	LoGiCJawY1_3,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);


  G4double distance_Source_JawY1_4=248.5*mm+JawY1_4_Dim_z;
  G4double dy_y1_4=distance_Source_JawY1_4*tan(theta);
  G4RotationMatrix RMY1_4,RMY2_4;   // unit rotation matrix
  G4double angleY1_4 = -theta;   // rotational angle
  RMY1_4.rotateX(angleY1_4);
  RMY2_4.rotateX(-angleY1_4);
  G4double jawY1_4_centre_y=JawY1_4_Dim_y+dy_y1_4;
  G4double jawY1_4_centre_x=0*mm;
 G4double jawY1_4_centre_z=248.5*mm+JawY1_4_Dim_z+dec;
 // G4double jawY1_4_centre_z=jawY1_3_centre_z-JawY1_3_Dim_z-JawY1_4_Dim_z;
  pLoGiCJawY1_4= new G4PVPlacement(G4Transform3D(RMY1_4,	//rotation 
		 G4ThreeVector(jawY1_4_centre_x, jawY1_4_centre_y, jawY1_4_centre_z)),
		 "pLoGiCJawY1_4",   //its name  (2nd constructor)
	LoGiCJawY1_4,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);
  pLoGiCJawY2_4= new G4PVPlacement(G4Transform3D(RMY2_4,	//rotation 
		 G4ThreeVector(jawY1_4_centre_x, -jawY1_4_centre_y, jawY1_4_centre_z)),
		 "pLoGiCJawY2_4",   //its name  (2nd constructor)
	LoGiCJawY1_4,         //its logical volume 
		 this->pMere,              //its mother volume 
		 false,                 //no voidean operation 
		 0);


}

