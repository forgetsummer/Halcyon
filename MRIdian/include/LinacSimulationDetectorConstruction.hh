#ifndef LinacSimulationDetectorConstruction_h
#define LinacSimulationDetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"
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
#include "projectDataTypeGlobals.hh"
#include "CellLayoutInitializer.hh"

#include "DicomDetectorConstruction.hh"

class LinacSimulationDetectorConstruction:public G4VUserDetectorConstruction

{
  
  public:
	 LinacSimulationDetectorConstruction(); // this is the constructor the class
	~LinacSimulationDetectorConstruction(); // this is the destructor of the class

        G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	G4ThreeVector GetWaterTankDimension()const;
	G4ThreeVector GetWaterTankCenterPosition()const;
	G4ThreeVector GetDoseTallyVolumeCenterPosition(G4String volumeType)const;
	G4ThreeVector GetCellDimension(G4String volumeType)const;
	CellLayoutInitializer GetCellLayoutInfo() const;
	
	std::vector<double>const& GetCellPositionXVec(G4String volumeType) const;
	std::vector<double>const& GetCellPositionYVec(G4String volumeType)const;
	std::vector<double>const& GetCellPositionZVec(G4String volumeType)const;
	
	G4ThreeVector GetDoseTallyVolumeDimension(G4String volumeType) const;
	
	G4ThreeVector GetVoxelNumberXYZ() const;
	G4ThreeVector GetVoxelDimension() const;
	G4ThreeVector GetVoxelDimMinXYZ() const;
	G4ThreeVector GetVoxelDimMaxXYZ() const;
	

private:

	void WaterTank(G4double Z_topSurface);
	void WaterTankHalcyon();
	void DoseTallyVolume(G4String volumeType,  G4ThreeVector centerPos, G4double dimX, double dimY, double dimZ);
	void DoseTallyVolumeVoxlization(G4String volumeType, G4int numX, G4int numY, G4int numZ);
	void SetDepthForDoseProfileXYTally(G4double depth);
	G4double waterTankDimX;
	G4double waterTankDimY;
	G4double waterTankDimZ;
	G4double waterTankPos_x;
	G4double waterTankPos_y;
	G4double waterTankPos_z;
	G4double cellHomeX;
	G4double cellHomeY;
	G4double cellHomeZ;
	
	G4double cellHomeX1;
	G4double cellHomeY1;
	G4double cellHomeZ1;
	
	G4double doseTallyVolumeDimX;
	G4double doseTallyVolumeDimY;
	G4double doseTallyVolumeDimZ;
	
	G4double doseTallyVolumeDimX1;
	G4double doseTallyVolumeDimY1;
	G4double doseTallyVolumeDimZ1;

	
	G4ThreeVector centerOfDoseTallyVolume;
	G4ThreeVector centerOfDoseTallyVolume1;
	
	G4VPhysicalVolume* waterTank_phys;
	G4VPhysicalVolume* PVWorld;
	G4Region* waterTankPhantom;
	G4Region* DICOMPhantomRegion;
	CellLayoutInitializer cellLayout;
	std::vector<double>cellPositionX; // define a vector to store the cell position
	std::vector<double>cellPositionY;
	std::vector<double>cellPositionZ;
	
	// later on, we may obtain the pdd and profile simultaneously, so we have these coordinates vectors
	CellLayoutInitializer cellLayout1;
	std::vector<double>cellPositionX1; // define a vector to store the cell position
	std::vector<double>cellPositionY1;
	std::vector<double>cellPositionZ1;
	
	G4LogicalVolume* BuildPhantomFromDICOM(); // function for building phantom from DICOm images
	void DICOMPhantom();
	
	DicomDetectorConstruction* theGeometry = 0;
};



#endif