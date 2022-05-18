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
// $Id$
//
/// \file B3aEventAction.cc
/// \brief Implementation of the B3aEventAction class

#include "B3aEventAction.hh"
#include "B3aRunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunCollectTriggeredSignals.hh"
#include "CellLayoutInitializer.hh"
#include "LinacSimulationDetectorConstruction.hh"
#include "G4IAEAphspWriter.hh"

#include "ArgumentInterpreter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction(B3aRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_cryst(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event* evt)
{
//   G4cout<<"get here"<<G4endl;
//  G4cout<<"The processing event now is: "<<evt->GetEventID()<<G4endl;
  currentEventID = evt->GetEventID();
   if (evt->GetEventID()%100000 ==0)
  {
    G4cout<<"The processing event now is: "<<evt->GetEventID()<<G4endl;
  }
  
  if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "write")
  {
     G4IAEAphspWriter* IAEAWriter = G4IAEAphspWriter::GetInstance();
    IAEAWriter->BeginOfEventAction(evt); // calling the function of G4IAEAphspWriter
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event* evt )
{
//   cerr<<"get to EndOfEventAction"<<endl;
//   PrintEDEPInformation(evt->GetEventID());
  

  MappingEDEPInfoToDoseVector(); // use this to tally the dose of water phantom
  MappingDicomDoseInfoToDoseVector(evt); // use this to tally the dose of dicom phantom
  
  
//   PrintCellDoseInformation(evt->GetEventID());
  
      
    std::map<G4int,G4double>::iterator mitr_cell;
    for (mitr_cell = cellEdepMap.begin();mitr_cell!=cellEdepMap.end();mitr_cell++)
    {
       G4RunCollectTriggeredSignals::CollectTriggeredSignal(mitr_cell->first,mitr_cell->second); // collecting total deposited energy here
    }
    
    
    for (mitr_cell = cellEdepMap1.begin(); mitr_cell != cellEdepMap1.end(); mitr_cell ++)

    {
      G4RunCollectTriggeredSignals::CollectTriggeredSignalPrimary(mitr_cell->first, mitr_cell->second);
    }
    
    for (mitr_cell = cellEdepMap2.begin(); mitr_cell!=cellEdepMap2.end(); mitr_cell ++)
    {
      G4RunCollectTriggeredSignals::CollectTriggeredSignalScatter(mitr_cell->first,mitr_cell->second);
    }
    
  edepVector.clear();
  cellEdepMap.clear();
  cellEdepMap1.clear();
  cellEdepMap2.clear();

 

}  


void B3aEventAction::PrintEDEPInformation(int eventID)
{

     if (edepVector.size()>0)
     {
       cout<<"The edep information of event: "<<eventID<<" is as "<<endl;
    }

  for (int i = 0 ;i<edepVector.size(); i++)
  {
    cout<<"x is: "<<edepVector[i].x<<" y is: "<<edepVector[i].y<<" z is: "<<edepVector[i].z<< " edep is: "<<edepVector[i].edep<<endl;
  }
  

}

void B3aEventAction::CollectEDEPInformation(EDEPInformation edepInfo)
{
  edepVector.push_back(edepInfo);

}

void B3aEventAction::MappingEDEPInfoToDoseVector()
{

    const LinacSimulationDetectorConstruction* detector= static_cast< const LinacSimulationDetectorConstruction*> \
  (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  
  G4String volumeType = "vertical";

  const std::vector<double>& cellPositionXVec = detector->GetCellPositionXVec(volumeType);
  const std::vector<double>& cellPositionYVec = detector->GetCellPositionYVec(volumeType);
  const std::vector<double>& cellPositionZVec = detector->GetCellPositionZVec(volumeType);
  
  
 
  G4ThreeVector cellHomeDimension = detector->GetCellDimension(volumeType);
  G4double cellHomeX = cellHomeDimension.getX();
  G4double cellHomeY = cellHomeDimension.getY();
  G4double cellHomeZ = cellHomeDimension.getZ();
  
  
  G4ThreeVector doseTallVolumeDimension = detector->GetDoseTallyVolumeDimension(volumeType);
  G4double doseTallyVolumeDimX = doseTallVolumeDimension.getX();
  G4double doseTallyVolumeDimY = doseTallVolumeDimension.getY();
  G4double doseTallyVolumeDimZ = doseTallVolumeDimension.getZ();
  
  G4int N_X = doseTallyVolumeDimX/cellHomeX;
  G4int N_Y = doseTallyVolumeDimY/cellHomeY;
  G4int N_Z = doseTallyVolumeDimZ/cellHomeZ;
  
  G4double doseTallVolumeISO_X = detector->GetDoseTallyVolumeCenterPosition(volumeType).getX();
  G4double doseTallVolumeISO_Y = detector->GetDoseTallyVolumeCenterPosition(volumeType).getY();
  G4double doseTallVolumeISO_Z = detector->GetDoseTallyVolumeCenterPosition(volumeType).getZ();
//   
  
  
  volumeType = "lateral";
  
  const std::vector<double>& cellPositionXVec1 = detector->GetCellPositionXVec(volumeType);
  const std::vector<double>& cellPositionYVec1 = detector->GetCellPositionYVec(volumeType);
  const std::vector<double>& cellPositionZVec1 = detector->GetCellPositionZVec(volumeType);
  

  G4ThreeVector cellHomeDimension1 = detector->GetCellDimension(volumeType);
  G4double cellHomeX1 = cellHomeDimension1.getX();
  G4double cellHomeY1 = cellHomeDimension1.getY();
  G4double cellHomeZ1 = cellHomeDimension1.getZ();
  
  G4ThreeVector doseTallVolumeDimension1 = detector->GetDoseTallyVolumeDimension(volumeType);
  G4double doseTallyVolumeDimX1 = doseTallVolumeDimension1.getX();
  G4double doseTallyVolumeDimY1 = doseTallVolumeDimension1.getY();
  G4double doseTallyVolumeDimZ1 = doseTallVolumeDimension1.getZ();
  
  G4int N_X1 = doseTallyVolumeDimX1/cellHomeX1;
  G4int N_Y1 = doseTallyVolumeDimY1/cellHomeY1;
  G4int N_Z1 = doseTallyVolumeDimZ1/cellHomeZ1;
  
  
  
  
  G4double doseTallVolumeISO_X1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getX();
  G4double doseTallVolumeISO_Y1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getY();
  G4double doseTallVolumeISO_Z1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getZ();
  
 
  
//   cout<<"doseTallVolumeISO_X = "<<doseTallVolumeISO_X<<endl;
//   cout<<"doseTallVolumeISO_Y = "<<doseTallVolumeISO_Y<<endl;
//   cout<<"doseTallVolumeISO_Z = "<<doseTallVolumeISO_Z<<endl;
//   
//   cout<<"cellHomeX = "<<cellHomeX<<" cellHomeY = "<<cellHomeY<<" cellHomeZ = "<<cellHomeZ<<endl;
//   for (int i = 0; i<cellPositionZVec.size(); i++)
//   {
//     cout<<"cell ID: "<<i<< " z= "<<cellPositionZVec[i] + doseTallVolumeISO_Z
//     <<" x= "<<cellPositionXVec[i] + doseTallVolumeISO_X<<" y= "<<cellPositionYVec[i] + doseTallVolumeISO_Y <<endl;
//   }
  
//   cout<<"doseTallVolumeISO_X1 = "<<doseTallVolumeISO_X1<<endl;
//   cout<<"doseTallVolumeISO_Y1 = "<<doseTallVolumeISO_Y1<<endl;
//   cout<<"doseTallVolumeISO_Z1 = "<<doseTallVolumeISO_Z1<<endl;
//   
//   cout<<"cellHomeX1 = "<<cellHomeX1<<" cellHomeY1 = "<<cellHomeY1<<" cellHomeZ1 = "<<cellHomeZ1<<endl;
//   for (int i = 0; i<cellPositionZVec1.size(); i++)
//   {
//     cout<<"cell ID: "<<i<< " z= "<<cellPositionZVec1[i] <<" x= "<<cellPositionXVec1[i]<<" y= "<<cellPositionYVec1[i] <<endl;
//   }
//  
//   
//   cout<<"N_X = "<<N_X<<" N_Y = "<<N_Y<<" N_Z = "<<N_Z<<endl;
//   cout<<"N_X1 = "<<N_X1<<" N_Y1 = "<<N_Y1<<" N_Z1 = "<<N_Z1<<endl;
  
//   exit(0);
  
  
  
  
  //////////////////////////////////////////////////////////////
  ////////////////for testing the phase space file distribution
//     for (int i = 0; i<edepVector.size();i++) // loop the edep vector of current event
//   {
//     for (int cellID = 0; cellID<cellPositionXVec.size(); cellID ++ )
//     {
// //       G4double x_cell = cellPositionXVec[cellID] + waterTankISO_X;
// //       G4double y_cell = cellPositionYVec[cellID] + waterTankISO_Y;
// //       G4double z_cell = cellPositionZVec[cellID] + waterTankISO_Z;
//       G4double x_cell = cellPositionXVec[cellID] + doseTallVolumeISO_X;
//       G4double y_cell = cellPositionYVec[cellID] + doseTallVolumeISO_Y;
//       G4double z_cell = cellPositionZVec[cellID] + doseTallVolumeISO_Z;
// //       cout<<"z_cell "<<z_cell<<endl;
//       G4double x = edepVector[i].x;
//       G4double y = edepVector[i].y;
//       G4double z = edepVector[i].z;
//       G4double edep = edepVector[i].edep;
//       if ((x_cell-cellHomeX/2< x && x <= x_cell+ cellHomeX/2) && (y_cell-cellHomeY/2 <y && y<=y_cell+cellHomeY/2)) // if this edep point is inside i-th voxel
//       {
// // 	cout<<"Inside cell "<<cellID<< " cell coordinate: "<<"x="<<x_cell<<" "<<"y="<<y_cell<<" "<<"z="<<z_cell<<" cellHomeX="<<cellHomeX
// // 	<<" cellHomeY="<<cellHomeY<<" cellHomeZ="<<cellHomeZ<<endl;
// 	cellEdepMap[cellID] = cellEdepMap[cellID] + 1;
// // 	cout<<"z = "<<z<<endl;
// 	break;
//       }
// 
//     }
//   }
// 
//   
  
  
  
//   cout<<"The voxel size is "<<cellPositionXVec.size()<<endl;
//   exit(0);
  
  ///////////////////////////////////////////////////////////////////
  
    
  for (int i = 0; i<edepVector.size();i++) // loop the edep vector of current event
  {
    
      G4double x = edepVector[i].x - doseTallVolumeISO_X; // coordinate translation back to dose tally volume's coordinate system
      G4double y = edepVector[i].y - doseTallVolumeISO_Y;
      G4double z = edepVector[i].z - doseTallVolumeISO_Z;
      G4double edep = edepVector[i].edep;

//       cout<<"x + doseTallyVolumeDimX - cellHomeX/2.0)/cellHomeX "<<(x + doseTallyVolumeDimX - cellHomeX/2.0)/cellHomeX<<endl;
      
      bool inX = (x>=-doseTallyVolumeDimX/2.0 && x<=doseTallyVolumeDimX/2.0);
      bool inY = (y>=-doseTallyVolumeDimY/2.0 && y<=doseTallyVolumeDimY/2.0);
      bool inZ = (z>=-doseTallyVolumeDimZ/2.0 && z<=doseTallyVolumeDimZ/2.0);
      
      if (inX && inY && inZ)
      {
// 	cout<<"x = "<<x <<"  y = "<<y<<" z = "<<z <<endl;
// 	cout<<"x = "<<x + doseTallVolumeISO_X<<"  y = "<<y+doseTallVolumeISO_Y<<" z = "<<z+doseTallVolumeISO_Z <<endl;
// 	G4int ii  = floor((x + doseTallyVolumeDimX/2.0 - cellHomeX/2.0)/cellHomeX);
// 	G4int jj  = floor((y + doseTallyVolumeDimY/2.0 - cellHomeY/2.0)/cellHomeY);
// 	G4int kk  = floor((z + doseTallyVolumeDimZ/2.0 - cellHomeZ/2.0)/cellHomeZ);
	
	G4int ii  = floor((x + doseTallyVolumeDimX/2.0 )/cellHomeX);
	G4int jj  = floor((y + doseTallyVolumeDimY/2.0 )/cellHomeY);
	G4int kk  = floor((z + doseTallyVolumeDimZ/2.0 )/cellHomeZ);
	
//         cout<<"ii = "<<ii<<" jj = "<<jj<<" kk = "<<kk<<endl;
//         cout<<"cellID = "<<cellID<<endl;
	if (ii>=N_X) // handles the boudary condition
	{
	  ii = ii - 1;
	}
	if (jj>=N_Y)
	{
	  jj = jj -1;
	}
	if(kk>=N_Z)
	{
	  kk = kk - 1;
	}
	
	G4int cellID = N_X*N_Y*kk + N_Y*ii + jj;
	
// 	        cout<<"ii = "<<ii<<" jj = "<<jj<<" kk = "<<kk<<endl;
//         cout<<"cellID = "<<cellID<<endl;
// 	if (cellID>999)
// 	{
// 	  cout<<"x = "<<x <<"  y = "<<y<<" z = "<<z <<endl;
// 	  cout<<"x = "<<x + doseTallyVolumeDimX/2.0<<"  y = "<<y+ doseTallyVolumeDimY/2.0<<" z = "<<z + doseTallyVolumeDimZ/2.0<<endl;
// 	  cout<<"ii = "<<ii<<" jj = "<<jj<<" kk = "<<kk<<endl;
// 	  cout<<"cellID = "<<cellID<<endl;
// 	}
	
	cellEdepMap[cellID] = cellEdepMap[cellID] + edep; // if the edep point is inside this cell, then add the energy into this cell
      
      }
    
  }
  
//   for (int i = 0; i<edepVector.size();i++ )
//   {
//       G4double x = edepVector[i].x - doseTallVolumeISO_X1; // coordinate translation back to dose tally volume's coordinate system
//       G4double y = edepVector[i].y - doseTallVolumeISO_Y1;
//       G4double z = edepVector[i].z - doseTallVolumeISO_Z1;
//       G4double edep = edepVector[i].edep;
//       
//       bool inX = (x>=-doseTallyVolumeDimX1/2.0 && x<=doseTallyVolumeDimX1/2.0);
//       bool inY = (y>=-doseTallyVolumeDimY1/2.0 && y<=doseTallyVolumeDimY1/2.0);
//       bool inZ = (z>=-doseTallyVolumeDimZ1/2.0 && z<=doseTallyVolumeDimZ1/2.0);
//       
// 
//       if (inX && inY && inZ)
//       {
// // 	cout<<"doseTallyVolumeDimX1 "<<doseTallyVolumeDimX1 <<
// //       " doseTallyVolumeDimY1 "<<doseTallyVolumeDimY1<<
// //       " doseTallyVolumeDimZ1 "<<doseTallyVolumeDimZ1<<endl;
// // 	cout<<"x = "<<x <<"  y = "<<y<<" z = "<<z <<endl;
// // 	G4int ii1 = floor((x + doseTallyVolumeDimX1/2.0 - cellHomeX1/2.0)/cellHomeX1);
// // 	G4int jj1 = floor((y + doseTallyVolumeDimY1/2.0 - cellHomeY1/2.0)/cellHomeY1);
// // 	G4int kk1 = floor((z + doseTallyVolumeDimZ1/2.0 - cellHomeZ1/2.0)/cellHomeZ1);
// 	
// 	G4int ii1 = floor((x + doseTallyVolumeDimX1/2.0 )/cellHomeX1);
// 	G4int jj1 = floor((y + doseTallyVolumeDimY1/2.0 )/cellHomeY1);
// 	G4int kk1 = floor((z + doseTallyVolumeDimZ1/2.0 )/cellHomeZ1);
// 	
// 	if (ii1>=N_X1)// handles the boundary condition
// 	{
// 	  ii1 = ii1 -1;
// 	}
// 	if (jj1>N_Y1)
// 	{
// 	  jj1 = jj1-1;
// 	}
// 	if(kk1>N_Z1)
// 	{
// 	  kk1 = kk1 -1;
// 	}
// 	
// 	G4int cellID1 = N_X1*N_Y1*kk1 + N_Y1*ii1 + jj1;
// // 	cout<<"ii1 = "<<ii1<<" jj1 = "<<jj1<<" kk1 = "<<kk1<<endl;
// // 	cout<<"cellID1 = "<<cellID1<<endl;
// 	cellEdepMap1[cellID1] = cellEdepMap1[cellID1] + edep; // if the edep point is inside this cell, then add the energy into this cell
//       }
//   }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Firtt version design, loop all the voxels to find the index of voxel
  ///////////////////////// It turns that this design could be very slow when the number of voxels are big
  
//   for (int i = 0; i<edepVector.size();i++) // loop the edep vector of current event
//   {
//     for (int cellID = 0; cellID<cellPositionXVec.size(); cellID ++ )
//     {  
// //       cout<<"Finish voxel "<<cellID<<endl;
// //       G4double x_cell = cellPositionXVec[cellID] + waterTankISO_X;
// //       G4double y_cell = cellPositionYVec[cellID] + waterTankISO_Y;
// //       G4double z_cell = cellPositionZVec[cellID] + waterTankISO_Z;
//       G4double x_cell = cellPositionXVec[cellID] + doseTallVolumeISO_X;
//       G4double y_cell = cellPositionYVec[cellID] + doseTallVolumeISO_Y;
//       G4double z_cell = cellPositionZVec[cellID] + doseTallVolumeISO_Z;
// //       cout<<"z_cell "<<z_cell<<endl;
//       G4double x = edepVector[i].x;
//       G4double y = edepVector[i].y;
//       G4double z = edepVector[i].z;
//       G4double edep = edepVector[i].edep;
//       if ((x_cell-cellHomeX/2< x && x <= x_cell+ cellHomeX/2) && (y_cell-cellHomeY/2 <y && y<=y_cell+cellHomeY/2) && 
// 	(z_cell- cellHomeZ/2 < z && z <= z_cell+ cellHomeZ/2)) // if this edep point is inside i-th voxel
//       {
// // 	cout<<"Inside cell "<<cellID<< " cell coordinate: "<<"x="<<x_cell<<" "<<"y="<<y_cell<<" "<<"z="<<z_cell<<" cellHomeX="<<cellHomeX
// // 	<<" cellHomeY="<<cellHomeY<<" cellHomeZ="<<cellHomeZ<<endl;
// 	cellEdepMap[cellID] = cellEdepMap[cellID] + edep; // if the edep point is inside this cell, then add the energy into this cell
// // 	cout<<"z = "<<z<<endl;
// 	break;
//       }
// 
//     }
//   }
//   
//       
//   for (int i = 0; i<edepVector.size();i++) // loop the edep vector of current event
//   {
//     for (int cellID = 0; cellID<cellPositionXVec1.size(); cellID ++ )
//     {
// //       cout<<"Finish voxel "<<cellID<<endl;
// //       G4double x_cell = cellPositionXVec[cellID] + waterTankISO_X;
// //       G4double y_cell = cellPositionYVec[cellID] + waterTankISO_Y;
// //       G4double z_cell = cellPositionZVec[cellID] + waterTankISO_Z;
//       G4double x_cell = cellPositionXVec1[cellID] + doseTallVolumeISO_X1;
//       G4double y_cell = cellPositionYVec1[cellID] + doseTallVolumeISO_Y1;
//       G4double z_cell = cellPositionZVec1[cellID] + doseTallVolumeISO_Z1;
// //       cout<<"z_cell "<<z_cell<<endl;
//       G4double x = edepVector[i].x;
//       G4double y = edepVector[i].y;
//       G4double z = edepVector[i].z;
//       G4double edep = edepVector[i].edep;
//       if ((x_cell-cellHomeX1/2< x && x <= x_cell+ cellHomeX1/2) && (y_cell-cellHomeY1/2 <y && y<=y_cell+cellHomeY1/2) && 
// 	(z_cell- cellHomeZ1/2 < z && z <= z_cell+ cellHomeZ1/2)) // if this edep point is inside i-th voxel
//       {
// // 	cout<<"Inside cell "<<cellID<< " cell coordinate: "<<"x="<<x_cell<<" "<<"y="<<y_cell<<" "<<"z="<<z_cell<<" cellHomeX="<<cellHomeX
// // 	<<" cellHomeY="<<cellHomeY<<" cellHomeZ="<<cellHomeZ<<endl;
// 	cellEdepMap1[cellID] = cellEdepMap1[cellID] + edep; // if the edep point is inside this cell, then add the energy into this cell
// // 	cout<<"z = "<<z<<endl;
// 	break;
//       }
// 
//     }
//   }
  
}

void B3aEventAction::MappingDicomDoseInfoToDoseVector(const G4Event* evt )
{
     //Hits collections

//   cout<<"The processing event now is: "<<evt->GetEventID()<<endl;
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
               
   // Get hits collections IDs
  if (fCollID_cryst < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollID_cryst   = SDMan->GetCollectionID("phantomSD/DoseDeposit");
  }
  
  
  std::map<G4int,G4double*>::iterator itr; // since it is a evtMap is double value map, so here we declare this map iterator
  // it is worth noting that the key of this map is the cell ID, or the copy number of the logical volume in this simulation.
  // so the key of this map will be the crystal detector, so we can use this key to decide which crystal is hitted by
  // the irradiating particle
   
  G4THitsMap<G4double>* evtMap = 
                     (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst)); // get the hit collection map by collection id
		     // since the energy is a double number, so here, the G4THitsMap is a double value map
//   cout<<"the size of evtMap is "<<evtMap->GetSize()<<endl;
//   exit(0);
		     
/*  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
  {
    ///G4int copyNb  = (itr->first);
//     G4cout<<"The detector ID is "<<(itr->first)<<G4endl;
    G4double edep = *(itr->second);
     G4cout<<"EvenID is "<<evt->GetEventID()<<" The dose "<<" in cellID "<<(itr->first)<<" is "<<edep<<G4endl;
  } */	

//     const LinacSimulationDetectorConstruction* detector= static_cast< const LinacSimulationDetectorConstruction*> \
//   (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//   
//     G4int fNVoxelX = detector->GetVoxelNumberXYZ().getX();
//     G4int fNVoxelY = detector->GetVoxelNumberXYZ().getY();
//     G4int fNVoxelZ = detector->GetVoxelNumberXYZ().getZ();
//     
//     G4double fVoxelDimX = detector->GetVoxelDimension().getX();
//     G4double fVoxelDimY = detector->GetVoxelDimension().getY();
//     G4double fVoxelDimZ = detector->GetVoxelDimension().getZ();
    
//     cout<<"fVoxelDimX= "<<fVoxelDimX <<" fVoxelDimY= "<<fVoxelDimY<<" fVoxelDimZ "<<fVoxelDimZ<<endl;
//      exit(0);
    
//     G4double fVoxelDimMinX = detector->GetVoxelDimMinXYZ().getX();
//     G4double fVoxelDimMinY = detector->GetVoxelDimMinXYZ().getY();
//     G4double fVoxelDimMinZ = detector->GetVoxelDimMinXYZ().getZ();
//     
//     G4double fVoxelDimMaxX = detector->GetVoxelDimMaxXYZ().getX();
//     G4double fVoxelDimMaxY = detector->GetVoxelDimMaxXYZ().getY();
//     G4double fVoxelDimMaxZ = detector->GetVoxelDimMaxXYZ().getZ();
    
    

//   for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
//   {
//     G4int copyNb  = (itr->first);
//     G4double x = fVoxelDimX*(copyNb%fNVoxelX);
//     G4double y = fVoxelDimX*((copyNb/fNVoxelX)%fNVoxelY);
//     G4double z = fVoxelDimZ*(copyNb/(fNVoxelX*fNVoxelY));
//     
//   }
    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
      G4int copyNb  = (itr->first);
//       cout<<"The dose voxel id is "<<copyNb<<endl;
      G4double dose = *(itr->second);
      cellEdepMap2[copyNb] = dose; // save dicom voxel dose information to cellEdepMap2
      
//       if (dose<0)
//       {G4cout<<"dose is "<<dose<<G4endl;}
      
    }
		     
}


void B3aEventAction::PrintCellDoseInformation(int eventID)
{
    if (edepVector.size()>0)
    {

	 cout<<"The cell edep information: "<<eventID<<" is as "<<endl;
   
    }
    
    std::map<G4int,G4double>::iterator mitr_cell;
    for (mitr_cell = cellEdepMap.begin();mitr_cell!=cellEdepMap.end();mitr_cell++)
    {
      cout<<"The cell is: "<<mitr_cell->first<<" Dose is: "<<mitr_cell->second<<endl;
    }
  

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
