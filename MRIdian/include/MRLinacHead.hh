#ifndef MRLinacHead_h
#define MRLinacHead_h 1

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ProductionCuts.hh"
#include "G4VUserDetectorConstruction.hh"
#include "projectDataTypeGlobals.hh"
#include "CellLayoutInitializer.hh"

// class CML2Acc1Messenger;
class MRLinacHead
{
public:
	 MRLinacHead(); // this is the constructor the class
	~MRLinacHead(); // this is the destructor of the class
	void GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume);
	G4ThreeVector GetTargetPosition() {return targetPosition;};

	

private:

	G4double maxConeAngle;
	G4double diameterOfPrimaryCollimatorTop;
	
	void Target();
	void PrimaryCollimator();
	void ElectronAbsorber();
	void IonizationChamber();
	void MLC();
	void MLCSteepEdge();
	void MLCCustomized(G4double coneAngle);
	void WaterTank();
	void WaterTankVoxlization(G4int numX,G4int numY,G4int numZ);
	void PhaseSpacePlane();
	G4double waterTankDimX;
	G4double waterTankDimY;
	G4double waterTankDimZ;
	G4double waterTankPos_x;
	G4double waterTankPos_y;
	G4double waterTankPos_z;
	G4double cellHomeX;
	G4double cellHomeY;
	G4double cellHomeZ;
	G4LogicalVolume *accWorldLV;
	G4VPhysicalVolume* PVWorld;
        G4VPhysicalVolume* targetTungsten_phys;
        G4VPhysicalVolume* targetCu_phys;
	G4VPhysicalVolume* primaryCollimator_phys;
	G4VPhysicalVolume* electronAbsorber_phys;
	G4VPhysicalVolume* ionizationChamber_phys;
	G4VPhysicalVolume* firstMLC_phys;
	G4VPhysicalVolume* secondMLC_phys;
	G4VPhysicalVolume* thirdMLC_phys;
	G4VPhysicalVolume* fourthMLC_phys;
	
	G4ThreeVector targetPosition;

};



#endif
