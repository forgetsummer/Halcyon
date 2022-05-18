#include "LinacSimulationDetectorConstruction.hh"
#include "MRLinacHead.hh"
#include "VarianClinac2100Head.hh"
#include "G4UserLimits.hh"
#include "TruebeamHead.hh"
#include "ArgumentInterpreter.hh"
#include "ReadLinacRTInformation.hh"
#include "G4GeometryManager.hh"

LinacSimulationDetectorConstruction::LinacSimulationDetectorConstruction()
{
  
}
LinacSimulationDetectorConstruction::~LinacSimulationDetectorConstruction()
{

}

G4VPhysicalVolume* LinacSimulationDetectorConstruction::Construct()
{
        // create the accelerator-world box
      G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      G4ThreeVector halfSize;

      G4double LinacWorldDimX = 2000*mm;
      G4double LinacWorldDimY = 2000*mm;
      G4double LinacWorldDimZ =2000*mm;
      halfSize.set(LinacWorldDimX,LinacWorldDimY,LinacWorldDimZ);
      
      
      G4GeometryManager::GetInstance()->SetWorldMaximumExtent(LinacWorldDimX*2); // this is for setting the geomety tolerance
      
      
      G4Box *accWorldB = new G4Box("accWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
      G4LogicalVolume *accWorldLV = new G4LogicalVolume(accWorldB, Vacuum, "accWorldL", 0, 0, 0);

     PVWorld =  new G4PVPlacement(0,                       //no rotation
	    G4ThreeVector(0,0,0),         //at (0,0,0)
	    accWorldLV,            //its logical volume
	    "linacWorldLV",               //its name
	    0,              //its mother  volume
	    false,                   //no boolean operation
	    0,                       //copy number
	    0);         // checking overlap*/
     
      G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
      simpleAlSVisAtt->SetVisibility(false);
// 		simpleAlSVisAtt->SetForceWireframe(false);
//       accWorldLV->SetVisAttributes(simpleAlSVisAtt);

      G4double z_target;
      
//       G4String LinacType = "MRLinac";
//       G4String LinacType = "VarianClinac2100";
      G4String LinacType = "Trubeam";
      
      if (LinacType == "MRLinac") // this is for MR Linac simulation
      {
	MRLinacHead* mrlinac = new MRLinacHead();
	mrlinac->GetLinacHeadLogicVolume(PVWorld);
	z_target = mrlinac->GetTargetPosition().getZ();
      }
      if (LinacType == "VarianClinac2100") // this is for VarianClinac2100  simulation 
      {     
	VarianClinac2100Head* varianLinac = new VarianClinac2100Head();
	varianLinac->GetLinacHeadLogicVolume(PVWorld);
	z_target = varianLinac->GetTargetPosition().getZ();// this is the z position of the tungsten target for generating photon
	
      }
      if (LinacType == "Trubeam")
      {
	TruebeamHead* truebeam = new TruebeamHead();
	truebeam->GetLinacHeadLogicVolume(PVWorld);
	z_target = truebeam->GetTargetPosition().getZ();
	
      }
      
      cout<<"z_target = "<<z_target<<endl;
//       exit(0);
      G4double SSD = 900*mm; // here we set up the SSD
      G4double z_waterTankTopSurface = z_target - SSD;
      cout<<"z_waterTankTopSurface= "<<z_waterTankTopSurface<<endl;
//       exit(0);
//       WaterTank(z_waterTankTopSurface);
//       WaterTankHalcyon(); // this is for calculating the water phantom dose
        DICOMPhantom(); // this is for calculating the dicom phantom dose
      
      
      G4ThreeVector doseTallyVolumePos = GetWaterTankCenterPosition();
      
      
      bool volumeVertical = true;
      bool volumeLateral =  false;
      if (volumeLateral)
      {
	G4String volumeType = "lateral";
	DoseTallyVolume(volumeType, doseTallyVolumePos,300*mm,5*mm,5*mm); // for profile simulation
	DoseTallyVolumeVoxlization(volumeType, 100,1,1);
	SetDepthForDoseProfileXYTally(100); //in unit of mm,  need this when we get profile at certain depth
      }
      

      if (volumeVertical)
      {
	G4String volumeType = "vertical";
	DoseTallyVolume(volumeType, doseTallyVolumePos,300*mm,300*mm,300*mm); // for geeting pdd
	DoseTallyVolumeVoxlization(volumeType, 120,120,120);
      }


  return PVWorld;

}

void LinacSimulationDetectorConstruction::ConstructSDandField()
{
  if (theGeometry) //we just use the sensitive detectors in dicom phantom volume
  {
   theGeometry->ConstructSDandField(); // here add the sensitive detector in DicomDetectorConstruction to LinacSimulationDetectorConstruction
  }
  
}


G4ThreeVector LinacSimulationDetectorConstruction::GetCellDimension(G4String volumeType) const
{
  if (volumeType == "vertical")
  {
    return G4ThreeVector(cellHomeX,cellHomeY,cellHomeZ);
  }
  if (volumeType == "lateral")
  {
    return G4ThreeVector(cellHomeX1,cellHomeY1,cellHomeZ1);
  }
}
CellLayoutInitializer LinacSimulationDetectorConstruction::GetCellLayoutInfo() const
{
    return cellLayout;

}

const vector< double >& LinacSimulationDetectorConstruction::GetCellPositionXVec(G4String volumeType) const
{
  

  if (volumeType == "vertical")
  {
    return cellPositionX;
  }
  if (volumeType == "lateral")
  {
    return cellPositionX1;
  }
  
}
const vector< double >& LinacSimulationDetectorConstruction::GetCellPositionYVec(G4String volumeType) const
{
  if (volumeType == "vertical")
  {
     return cellPositionY;
  }
  if (volumeType == "lateral")
  {
    return cellPositionY1;
  }

}

const vector< double >& LinacSimulationDetectorConstruction::GetCellPositionZVec(G4String volumeType) const
{
  if (volumeType == "vertical")
  {
    return cellPositionZ;
  }
  if (volumeType == "lateral")
  {
    return cellPositionZ1;
  }
  
}


G4ThreeVector LinacSimulationDetectorConstruction::GetWaterTankCenterPosition() const
{
  return G4ThreeVector(waterTankPos_x,waterTankPos_y,waterTankPos_z);
}
G4ThreeVector LinacSimulationDetectorConstruction::GetWaterTankDimension() const
{
   return G4ThreeVector(waterTankDimX,waterTankDimY,waterTankDimZ);
}

void LinacSimulationDetectorConstruction::WaterTank(G4double Z_topSurface)
{
      G4Material *water=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
      
//       G4Material *water=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"); // for checking the spatial distribution of phase space file
      waterTankDimX = 400*mm;
      waterTankDimY = 400*mm;
      waterTankDimZ = 400*mm;
      G4Box* waterTank = new G4Box("waterTank",waterTankDimX/2,waterTankDimY/2,waterTankDimZ/2);
      G4LogicalVolume* waterTankLogicV = new G4LogicalVolume(waterTank,water,"waterTankLogicV",0,0,0);
      
      waterTankPos_x = 0.*mm;
      waterTankPos_y = 0.*mm;
      waterTankPos_z = Z_topSurface-waterTankDimZ/2;
      
      waterTank_phys = new G4PVPlacement(0,G4ThreeVector(waterTankPos_x,waterTankPos_y,waterTankPos_z),"waterTank_phys",waterTankLogicV,
	PVWorld,false,0
      );

      G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
      simpleAlSVisAtt->SetVisibility(true);
      waterTankLogicV->SetVisAttributes(simpleAlSVisAtt);
      
//       waterTankPhantom->AddRootLogicalVolume(waterTankLogicV);
      
//       G4Region *regVol;
//       regVol= new G4Region("targetR");
//       G4ProductionCuts* cuts = new G4ProductionCuts;
//       cuts->SetProductionCut(0.05*mm); // set up production cuts, mainly for secondary particles
//       regVol->SetProductionCuts(cuts);

      
//       waterTankLogicV->SetUserLimits(new G4UserLimits(0.005*mm));
      
}

void LinacSimulationDetectorConstruction::WaterTankHalcyon()
{
        G4Material *water=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
      
//       G4Material *water=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"); // for checking the spatial distribution of phase space file
      waterTankDimX = 300*mm;
      waterTankDimY = 300*mm;
      waterTankDimZ = 300*mm;
      G4Box* waterTank = new G4Box("waterTank",waterTankDimX/2,waterTankDimY/2,waterTankDimZ/2);
      G4LogicalVolume* waterTankLogicV = new G4LogicalVolume(waterTank,water,"waterTankLogicV",0,0,0);
      
     G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
     ReadLinacRTInformation readHlcy;
     readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information
     
     G4double couchLat = readHlcy.GetCoutchInformation().getX()*cm;
     G4double couchVrt = readHlcy.GetCoutchInformation().getZ()*cm;
     G4double couchLng = readHlcy.GetCoutchInformation().getY()*cm;
     
     G4double waterTankDimXTPS = 30*cm;
     G4double waterTankDimYTPS = 30*cm;
     G4double waterTankDimZTPS = 30*cm;
     
     waterTankPos_x = 1/2.0*waterTankDimXTPS - 1/2.0*waterTankDimX + couchLat;
     waterTankPos_y = 1/2.0*waterTankDimYTPS - 1/2.0*waterTankDimY + couchLng;
     waterTankPos_z = 1/2.0*waterTankDimZTPS - 1/2.0*waterTankDimZ + couchVrt;
     cout<<"1/2waterTankDimX "<<1/2*waterTankDimX<<endl;
     
     cout<<"CouchLat: "<<couchLat<< " CouchVrt: "<<couchVrt<<" CouchLng: "<<couchLng<<endl;
     cout<<"waterTankPos_x= "<<waterTankPos_x<<" waterTankPos_y= "<<waterTankPos_y<<
     " waterTankPos_z= "<<waterTankPos_z<<endl;
//      exit(0);
     
 
//       waterTankPos_x = 0.*mm;
//       waterTankPos_y = 0.*mm;
//       waterTankPos_z = Z_topSurface-waterTankDimZ/2;
      
      waterTank_phys = new G4PVPlacement(0,G4ThreeVector(waterTankPos_x,waterTankPos_y,waterTankPos_z),"waterTank_phys",waterTankLogicV,
	PVWorld,false,0
      );

      G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
      simpleAlSVisAtt->SetVisibility(true);
      waterTankLogicV->SetVisAttributes(simpleAlSVisAtt);
      
//       waterTankPhantom->AddRootLogicalVolume(waterTankLogicV);
      
//       G4Region *regVol;
//       regVol= new G4Region("targetR");
//       G4ProductionCuts* cuts = new G4ProductionCuts;
//       cuts->SetProductionCut(0.05*mm); // set up production cuts, mainly for secondary particles
//       regVol->SetProductionCuts(cuts);

      
//       waterTankLogicV->SetUserLimits(new G4UserLimits(0.005*mm));
      

}

void LinacSimulationDetectorConstruction::DICOMPhantom()
{
  G4LogicalVolume* logicalPatientFromDICOM = BuildPhantomFromDICOM();
  
  ///WARNING we tentatively use the watertank position to put the dicom phantom. Later we need to fin the 
  // right way to do it
  
      waterTankDimX = 300*mm;
      waterTankDimY = 300*mm;
      waterTankDimZ = 300*mm;
      
     G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
     ReadLinacRTInformation readHlcy;
     readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information
     
     G4double couchLat = readHlcy.GetCoutchInformation().getX()*cm;
     G4double couchVrt = readHlcy.GetCoutchInformation().getZ()*cm;
     G4double couchLng = readHlcy.GetCoutchInformation().getY()*cm;
     
     G4double waterTankDimXTPS = 30*cm;
     G4double waterTankDimYTPS = 30*cm;
     G4double waterTankDimZTPS = 30*cm;
     
     waterTankPos_x = 1/2.0*waterTankDimXTPS - 1/2.0*waterTankDimX + couchLat;
     waterTankPos_y = 1/2.0*waterTankDimYTPS - 1/2.0*waterTankDimY + couchLng;
     waterTankPos_z = 1/2.0*waterTankDimZTPS - 1/2.0*waterTankDimZ + couchVrt;
     cout<<"1/2waterTankDimX "<<1/2.0*waterTankDimX<<endl;
     
     cout<<"CouchLat: "<<couchLat<< " CouchVrt: "<<couchVrt<<" CouchLng: "<<couchLng<<endl;
     
     
     cout<<"dicom phantom dimension in X is: "<<theGeometry->GetVoxelDimMaxXYZ().getX()-theGeometry->GetVoxelDimMinXYZ().getX()<<endl;
     cout<<"dicom phantom dimension in Y is: "<<theGeometry->GetVoxelDimMaxXYZ().getY()-theGeometry->GetVoxelDimMinXYZ().getY()<<endl;
     cout<<"dicom phantom dimension in Z is: "<<theGeometry->GetVoxelDimMaxXYZ().getZ()-theGeometry->GetVoxelDimMinXYZ().getZ()<<endl;
     cout<<"dicom voxel size in X is: "<<theGeometry->GetVoxelDimension().getX()<<endl;
     cout<<"dicom voxel size in Y is: "<<theGeometry->GetVoxelDimension().getY()<<endl;
     cout<<"dicom voxel size in Z is: "<<theGeometry->GetVoxelDimension().getZ()<<endl;
     cout<<"number of voxels in X is: "<<theGeometry->GetVoxelNumberXYZ().getX()<<endl;
     cout<<"number of voxels in Y is: "<<theGeometry->GetVoxelNumberXYZ().getY()<<endl;
     cout<<"number of voxels in Z is: "<<theGeometry->GetVoxelNumberXYZ().getZ()<<endl;
     
//      exit(0);
     
//       couchLat = couchLat + 1e-3*mm;
//       couchLng = couchLng + 1e-3*mm;
//       couchVrt = couchVrt + 1e-3*mm;
      
      cout<<"couchLat = "<<couchLat<<" couchLng = "<< couchLng<<" couchVrt = "<<couchVrt<<endl;
//       exit(0);
      
     
             G4RotationMatrix rotmBtmp = G4RotationMatrix();
	G4double gantryRtn = 90*deg; // rotate 90 degree to make sure the dicom phantom in right direction
	G4double collRtn = 180*deg;
      rotmBtmp.rotateX(gantryRtn);
//       rotmBtmp.rotateZ(collRtn);
//   G4ThreeVector ISOPos = G4ThreeVector(couchLat,couchLng, couchVrt);
  G4ThreeVector ISOPos = G4ThreeVector(couchLat,couchLng, couchVrt);
//    G4ThreeVector ISOPos = G4ThreeVector(0,0, 0);
  G4Transform3D transform = G4Transform3D(rotmBtmp,ISOPos);
     
      waterTank_phys = new G4PVPlacement(transform,"dicomPhantom_phys",logicalPatientFromDICOM,
	PVWorld,false,true
      );
     
      
//     G4LogicalVolume* testLV = PVWorld->GetLogicalVolume();
//     cout<<"The number of daughter logical volumes are "<<testLV->GetNoDaughters()<<endl;
//     
//     for (int i = 0;i<testLV->GetNoDaughters();i++)
//     {
//       cout<<"The name of the daughter logical volume is "<<testLV->GetDaughter(i)->GetName()<<endl;
//     }
//     
//     for (int i = 0;i<testLV->GetNoDaughters();i++)// set the position of linac components according to the needs
//     {
//       G4ThreeVector position = testLV->GetDaughter(i)->GetTranslation();
//       cout<<"The initial position information of daughter volume are as: "<<position.getX()<<" "<<position.getY()<<" "<<position.getZ()<<endl;
// //       position.setX(position.getX()+ relativePositionOfLinacHead.getX());
// //       position.setY(position.getY()+ relativePositionOfLinacHead.getY());
// //       position.setZ(position.getZ()+ relativePositionOfLinacHead.getZ());
// //       testLV->GetDaughter(i)->SetTranslation(position);
//     }
//     exit(0);
      
            
//       DICOMPhantomRegion->AddRootLogicalVolume(logicalPatientFromDICOM);
//       
//       G4Region *regVol;
//       regVol= new G4Region("dicom");
//       G4ProductionCuts* cuts = new G4ProductionCuts;
//       cuts->SetProductionCut(0.05*mm); // set up production cuts, mainly for secondary particles
//       regVol->SetProductionCuts(cuts);
//       
//       
//       logicalPatientFromDICOM->SetUserLimits(new G4UserLimits(0.005*mm));
 
}



void LinacSimulationDetectorConstruction::DoseTallyVolume(G4String volumeType,G4ThreeVector centerPos, G4double dimX, double dimY, double dimZ)
{
  
  if (volumeType == "vertical")
  {
    centerOfDoseTallyVolume = centerPos;
    doseTallyVolumeDimX = dimX;
    doseTallyVolumeDimY = dimY;
    doseTallyVolumeDimZ = dimZ;
  }
  
  if (volumeType == "lateral")
  {
    centerOfDoseTallyVolume1 = centerPos;
    doseTallyVolumeDimX1 = dimX;
    doseTallyVolumeDimY1 = dimY;
    doseTallyVolumeDimZ1 = dimZ;
    
  }
  


}


void LinacSimulationDetectorConstruction::DoseTallyVolumeVoxlization(G4String volumeType,G4int numX, G4int numY, G4int numZ)
{
  if (volumeType == "vertical")
  {
    G4int doseTallyVolumeVoxlizeNumX = numX;
    G4int doseTallyVolumeVoxlizeNumY = numY;
    G4int doseTallyVolumeVoxlizeNumZ = numZ;
    cellHomeX = doseTallyVolumeDimX/doseTallyVolumeVoxlizeNumX;
    cellHomeY = doseTallyVolumeDimY/doseTallyVolumeVoxlizeNumY;
    cellHomeZ = doseTallyVolumeDimZ/doseTallyVolumeVoxlizeNumZ;
    
    
    
    cellLayout.SetCellHomeParamter(cellHomeX,cellHomeY,cellHomeZ);
    cellLayout.RectangularSlabFullMesh(doseTallyVolumeDimX,doseTallyVolumeDimY,doseTallyVolumeDimZ);
   
    
    cellPositionX = cellLayout.GetCellPositionXVec();
    cellPositionY = cellLayout.GetCellPositionYVec();
    cellPositionZ = cellLayout.GetCellPositionZVec();
    
    
    
//       for (int cellID = 0; cellID<cellLayout.GetCellNumber(); cellID ++ )
//       {
//         G4double x_cell = cellLayout.GetCellPositionX(cellID);
//         G4double y_cell = cellLayout.GetCellPositionY(cellID);
//         G4double z_cell = cellLayout.GetCellPositionZ(cellID);
//         
//         cout<<"cellID: "<<cellID<<" x= "<<x_cell<<" y= "<<y_cell<<" z= "<<z_cell<<endl;
//      
//       }
//       cout<<"cellHomeX = "<<cellHomeX<<" cellHomeY = "<<cellHomeY<<" cellHomeZ = "<<cellHomeZ<<endl;
//       cout<<"doseTallyVolumeDimX = "<<doseTallyVolumeDimX<<" doseTallyVolumeDimY = "<<doseTallyVolumeDimY 
//       <<" doseTallyVolumeDimZ = "<<doseTallyVolumeDimZ<<endl;
//       exit(0);
  }
  
  if (volumeType == "lateral")
  {
    G4int doseTallyVolumeVoxlizeNumX = numX;
    G4int doseTallyVolumeVoxlizeNumY = numY;
    G4int doseTallyVolumeVoxlizeNumZ = numZ;
    cellHomeX1 = doseTallyVolumeDimX1/doseTallyVolumeVoxlizeNumX;
    cellHomeY1 = doseTallyVolumeDimY1/doseTallyVolumeVoxlizeNumY;
    cellHomeZ1 = doseTallyVolumeDimZ1/doseTallyVolumeVoxlizeNumZ;
    
    
    
    cellLayout1.SetCellHomeParamter(cellHomeX1,cellHomeY1,cellHomeZ1);
    cellLayout1.RectangularSlabFullMesh(doseTallyVolumeDimX1,doseTallyVolumeDimY1,doseTallyVolumeDimZ1);
    
    
    cellPositionX1 = cellLayout1.GetCellPositionXVec();
    cellPositionY1 = cellLayout1.GetCellPositionYVec();
    cellPositionZ1 = cellLayout1.GetCellPositionZVec();
    
//       for (int cellID = 0; cellID<cellLayout.GetCellNumber(); cellID ++ )
//       {
//         G4double x_cell = cellLayout.GetCellPositionX(cellID);
//         G4double y_cell = cellLayout.GetCellPositionY(cellID);
//         G4double z_cell = cellLayout.GetCellPositionZ(cellID);
//         
//         cout<<"cellID: "<<cellID<<" x= "<<x_cell<<" y= "<<y_cell<<" z= "<<z_cell<<endl;
//      
//       }
//       cout<<"cellHomeX = "<<cellHomeX<<" cellHomeY = "<<cellHomeY<<" cellHomeZ = "<<cellHomeZ<<endl;
//       cout<<"doseTallyVolumeDimX = "<<doseTallyVolumeDimX<<" doseTallyVolumeDimY = "<<doseTallyVolumeDimY 
//       <<" doseTallyVolumeDimZ = "<<doseTallyVolumeDimZ<<endl;
//       exit(0);
  }



}



G4ThreeVector LinacSimulationDetectorConstruction::GetDoseTallyVolumeCenterPosition(G4String volumeType) const
{
  if (volumeType == "vertical")
  {
    return centerOfDoseTallyVolume;
  }
  if (volumeType == "lateral")
  {
    return centerOfDoseTallyVolume1;
  }
}


G4ThreeVector LinacSimulationDetectorConstruction::GetDoseTallyVolumeDimension(G4String volumeType) const
{
  if (volumeType == "vertical")
  {
    G4ThreeVector doseTallVolumeDimension(doseTallyVolumeDimX,doseTallyVolumeDimY, doseTallyVolumeDimZ);
    return doseTallVolumeDimension;
  }
  if(volumeType == "lateral")
  {
    G4ThreeVector doseTallVolumeDimension(doseTallyVolumeDimX1, doseTallyVolumeDimY1, doseTallyVolumeDimZ1);
    return doseTallVolumeDimension;
  }
  

}


void LinacSimulationDetectorConstruction::SetDepthForDoseProfileXYTally(G4double depth)
{
  for (int i = 0; i<cellPositionZ1.size(); i++)
  {
    cellPositionZ1[i] = cellPositionZ1[i] + waterTankDimZ/2.0 -  depth;
  }

}


////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
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
  
  

G4LogicalVolume* LinacSimulationDetectorConstruction::BuildPhantomFromDICOM()
{
  
  
  // here we start to modify the DICOM example to build the voxelized phantom
  
  G4bool bPartial = FALSE; // if we we need to use partial dicomo phatnom, then set this as true
  G4bool bNest = TRUE; // if we use the nested phantom, then set this as TRUE,
//   G4bool bNest = FALSE; // if we use the nested phantom, then set this as TRUE,
//    DicomDetectorConstruction* theGeometry = 0; // we put DicomDetectorConstruction as a private variable of LinacSimulationDetectorConstruction, so we comment this line
  
  G4String dicomFilePath = ArgumentInterpreter::GetDicomFilePath();
 
//   exit(0);
  
  #ifdef G4_DCMTK
    DicomFileMgr* theFileMgr = 0;
  #else
    DicomHandler* dcmHandler = 0;
  #endif
    
    if( !bPartial ) // if not use use the partial phantom
    {
  #ifdef G4_DCMTK // if we have DCMTK, then we do this
      
      theFileMgr = DicomFileMgr::GetInstance();
      theFileMgr->Convert("Data.dat");
      cout<<"we use DCMTK"<<endl;
      
  #else // else we do this, which is a regular way
      // Treatment of DICOM images before creating the G4runManager
//       cout<<"we use regular way to build geometry"<<endl;
//       exit(0);
      dcmHandler = new DicomHandler;
      
//       G4String dicomFilePath = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ProsIMRT/";
      cout<<"length of file is "<<sizeof(dicomFilePath)<<endl;
      cout<<"length in char is "<<sizeof(dicomFilePath)/sizeof(char);
//       exit(0);
//       dcmHandler->SetDICOMFilesFolderPath("/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ProsIMRT/");     
      dcmHandler->SetDICOMFilesFolderPath(dicomFilePath);     
      dcmHandler->CheckFileFormat();
  //     cout<<"Finishing dcmHandler specification"<<endl;
  //     exit(0);
  #endif
      // Initialisation of physics, geometry, primary particles ...
      if(bNest) // if use the nested phantom option
      {
	G4cout<<"Use DicomNestedParamDetectorConstruction"<<G4endl;
	theGeometry = new DicomNestedParamDetectorConstruction();
      }
      else // otherwise use the regular phantom option
      {
// 	cout<<"get here"<<endl;
// 	exit(0);
	G4cout<<"Use DicomRegularDetectorConstruction"<<G4endl;
	theGeometry = new DicomRegularDetectorConstruction();
      }
    } 
    else // otherwise build the partial dicom phantom
    {
      G4cout<<"Use DicomPartialDetectorConstruction"<<G4endl;
      theGeometry = new DicomPartialDetectorConstruction();
    }    
    
//     theGeometry->SetDICOMFilesFolderPath("/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ProsIMRT/");
       theGeometry->SetDICOMFilesFolderPath(dicomFilePath);
    cout<<"dicomFilePath is "<<dicomFilePath<<endl;
//      exit(0);
    theGeometry->Construct(); // Construct the geometry, inlcuding everything in DICOM example
    
    G4cout<<"Finishing building geometry"<<G4endl;

    G4cout<<"After initialize geometry in runManager"<<G4endl;
    
    
    G4VPhysicalVolume* thePhysicalVolume =  theGeometry->GetContainerPhysicalVolume();
    G4LogicalVolume* dicomLogicalVolume = theGeometry->GetContainerLogicalVolume();
//     G4LogicalVolume* dicomLogicalVolume = theGeometry->GetWorldLogicalVolume();
    
    
//     for (int i = 0;i<dicomLogicalVolume->GetNoDaughters();i++)
//     {
//        G4ThreeVector position =dicomLogicalVolume->GetDaughter(i)->GetTranslation();
//        cout<<"GetDaughter name is "<<dicomLogicalVolume->GetDaughter(i)->GetName()<<endl;
//        cout<<"position x "<<position.getX()<<","<<position.getY()<<","<<position.getZ()<<endl;
//        G4ThreeVector trans = dicomLogicalVolume->GetDaughter(i)->GetTranslation();
//        cout<<"position trans x "<<trans.getX()<<","<<trans.getY()<<","<<trans.getZ()<<endl;
//     }
//     exit(0);
    
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

G4ThreeVector LinacSimulationDetectorConstruction::GetVoxelNumberXYZ() const
{
  return theGeometry->GetVoxelNumberXYZ();
}
G4ThreeVector LinacSimulationDetectorConstruction::GetVoxelDimension() const
{
  return theGeometry->GetVoxelDimension();
}
G4ThreeVector LinacSimulationDetectorConstruction::GetVoxelDimMinXYZ() const
{
  return theGeometry->GetVoxelDimMinXYZ();
}
G4ThreeVector LinacSimulationDetectorConstruction::GetVoxelDimMaxXYZ() const
{
  return theGeometry->GetVoxelDimMaxXYZ();
}

  



