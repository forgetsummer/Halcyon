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
// $Id: B3aRunAction.cc 99559 2016-09-27 07:02:21Z gcosmo $
//
/// \file B3aRunAction.cc
/// \brief Implementation of the B3aRunAction class

#include "B3aRunAction.hh"
#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunCollectTriggeredSignals.hh"
#include "CellSimulationStatisticalErrorAnalysis.hh"
#include "ArgumentInterpreter.hh"
#include "LinacSimulationDetectorConstruction.hh"
#include "time.h"

#include "G4IAEAphspWriter.hh"

#include "ArgumentInterpreter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aRunAction::B3aRunAction()
 : G4UserRunAction(),
   fGoodEvents(0),
   fSumDose(0.)  
{  
  //add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fGoodEvents);
  accumulableManager->RegisterAccumulable(fSumDose); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aRunAction::~B3aRunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
     // IN RUNACTION CLASS:
    long seeds[2];
    long systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    G4Random::setTheSeeds(seeds);
  
 
  if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "write")
  {
      G4IAEAphspWriter* IAEAWriter = G4IAEAphspWriter::GetInstance();
      IAEAWriter->SetZStop(440*mm);  // here we set the position of phase space plane just above the MLC
      
      // we need to tell the G4IAEAphspWriter what is the primary particle direction
	
      IAEAWriter->SetPrimaryFlyDirection(G4ThreeVector(0,0,-1)); // here set up the primary photon flying direction
      
      IAEAWriter->SetFileName("MRLinacPhaseSpaceFile");
      G4IAEAphspWriter::GetInstance()->BeginOfRunAction(run);
  }

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B3PrimaryGeneratorAction* generatorAction
    = static_cast<const B3PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String partName;
  if (generatorAction) 
  {
    G4ParticleDefinition* particle 
      = generatorAction->GetParticleGun()->GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
          
  // Print results
  //
  if (IsMaster())
  {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------"
     << G4endl
     << "  The run was " << nofEvents << " events ";
  }
  else
  {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------"
     << G4endl
     << "  The run was " << nofEvents << " "<< partName;
  }      
  G4cout
     << "; Nb of 'good' e+ annihilations: " << fGoodEvents.GetValue()  << G4endl
     << " Total dose in patient : " << G4BestUnit(fSumDose.GetValue(),"Dose") 
     << G4endl 
     << "------------------------------------------------------------" << G4endl 
     << G4endl;
     
     if (IsMaster())
     {
       
	if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "write")
	{
	    G4IAEAphspWriter* IAEAWriter = G4IAEAphspWriter::GetInstance();
	    IAEAWriter->EndOfRunAction(run);
     
	}
  // here we do the final tally of the triggered signal of detectors
       
       
     
//           cout<<"the nofEvents is "<<nofEvents<<endl;
//           exit(0);
     
     G4RunCollectTriggeredSignals myCollection;
     std::map<int, G4double > cellEnergyMap = myCollection.GetWholeRunSignalCountMap();
     std::map<int, G4double > cellEnergySquaredMap = myCollection.GetWholeRunSignalCountSquaredMap();
     
     std::map<int, G4double> cellSignalMeanMap = GetMCResults(cellEnergyMap,cellEnergySquaredMap, nofEvents).first;
     std::map<int, G4double> cellSignalStdMap = GetMCResults(cellEnergyMap,cellEnergySquaredMap,nofEvents).second;
     
     
     // get  primary contribution
     std::map<int, G4double > cellEnergyPrimaryMap = myCollection.GetWholeRunSignalCountPrimaryMap();
     std::map<int, G4double > cellEnergyPrimarySquaredMap = myCollection.GetWholeRunSignalCountPrimarySquaredMap();
     
     
     std::map<int, G4double > cellSignalPrimaryMeanMap = GetMCResults(cellEnergyPrimaryMap,cellEnergyPrimarySquaredMap,nofEvents).first;
     std::map<int, G4double > cellSignalPrimaryStdMap = GetMCResults(cellEnergyPrimaryMap,cellEnergyPrimarySquaredMap,nofEvents).second;
     
     // get scatter contribution
     std::map<int, G4double > cellEnergyScatterMap = myCollection.GetWholeRunSignalCountScatterMap();
     std::map<int, G4double > cellEnergyScatterSquaredMap = myCollection.GetWholeRunSignalCountScatterSquaredMap();
     
     std::map<int, G4double > cellSignalScatterMeanMap = GetMCResults(cellEnergyScatterMap,cellEnergyScatterSquaredMap,nofEvents).first;
     std::map<int, G4double > cellSignalScatterStdMap  = GetMCResults(cellEnergyScatterMap,cellEnergyScatterSquaredMap,nofEvents).second;
     

     
     
//       G4String filename;
//       filename = ArgumentInterpreter::GetOutPutFileName();
      
      stringstream ss;
      ss<<run->GetRunID();
//       G4String outputDoseFileName = ss.str()+"_dose"+".csv";
      
//       G4String outputDoseFileName=filename+"_dose"+".csv";// write dose tally information out to a file
//       ofstream file_Dose;
      
//       file_Dose.open(outputDoseFileName);
      
      
          const LinacSimulationDetectorConstruction* detector= static_cast< const LinacSimulationDetectorConstruction*> \
      (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      
//       G4double waterTankISO_X = detector->GetWaterTankCenterPosition().getX();
//       G4double waterTankISO_Y = detector->GetWaterTankCenterPosition().getY();
//       G4double waterTankISO_Z = detector->GetWaterTankCenterPosition().getZ();
      
      G4String volumeType = "vertical";
      
      G4double doseTallVolumeISO_X = detector->GetDoseTallyVolumeCenterPosition(volumeType).getX();
      G4double doseTallVolumeISO_Y = detector->GetDoseTallyVolumeCenterPosition(volumeType).getY();
      G4double doseTallVolumeISO_Z = detector->GetDoseTallyVolumeCenterPosition(volumeType).getZ();
      
//       CellLayoutInitializer cellLayout = detector->GetCellLayoutInfo();
      const std::vector<double>& cellPositionXVec = detector->GetCellPositionXVec(volumeType);
      const std::vector<double>& cellPositionYVec = detector->GetCellPositionYVec(volumeType);
      const std::vector<double>& cellPositionZVec = detector->GetCellPositionZVec(volumeType);
      
      
      volumeType = "lateral";
      
      G4double doseTallVolumeISO_X1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getX();
      G4double doseTallVolumeISO_Y1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getY();
      G4double doseTallVolumeISO_Z1 = detector->GetDoseTallyVolumeCenterPosition(volumeType).getZ();
      
//       CellLayoutInitializer cellLayout = detector->GetCellLayoutInfo();
      const std::vector<double>& cellPositionXVec1 = detector->GetCellPositionXVec(volumeType);
      const std::vector<double>& cellPositionYVec1 = detector->GetCellPositionYVec(volumeType);
      const std::vector<double>& cellPositionZVec1 = detector->GetCellPositionZVec(volumeType);
      
     
      G4String filename = ArgumentInterpreter::GetOutPutFileName();
      G4String controlFileName = ArgumentInterpreter::GetControlPointsFile();
      controlFileName = controlFileName.substr(0,controlFileName.length()-5);
      
      
      G4String outputDoseFileName="PDD_"+ filename+ "_" +  controlFileName +  "_dose"+".csv";// write dose tally information out to a file
      
      G4String outputDoseFileName1 = "Profile_"+ filename + "_" +  controlFileName + "_dose" + ".csv";
      
      G4String outPutDoseFileName2 = "DICOM_" + filename + "_" +  controlFileName + "_dose" + ".csv";
      
      G4String voxelPhantomInfo = "VoxelPhantomInfo.csv";
      outputDoseFileName = ArgumentInterpreter::GetDoseFilePath() + outputDoseFileName;
      outputDoseFileName1 = ArgumentInterpreter::GetDoseFilePath() + outputDoseFileName1;
      outPutDoseFileName2 = ArgumentInterpreter::GetDoseFilePath() + outPutDoseFileName2;
      voxelPhantomInfo = ArgumentInterpreter::GetDoseFilePath() + voxelPhantomInfo;
      
      cout<<"outputDoseFileName "<<outputDoseFileName<<endl;
      cout<<"outputDoseFileName1 "<<outputDoseFileName1<<endl;
      cout<<"outPutDoseFileName2 "<<outPutDoseFileName2<<endl;
//        exit(0);
      G4String pathOfDoseFiles = ArgumentInterpreter::GetDoseFilePath();
      
      std::remove(outputDoseFileName); // first, remove the old file if the file exists
      std::remove(outputDoseFileName1);
      std::remove(outPutDoseFileName2);
      std::remove(voxelPhantomInfo);
      
      ofstream file_Dose;
      ofstream file_Dose1;
      ofstream file_Dose2;
      ofstream file_Dose3;
      
//       file_Dose.open(outputDoseFileName);
//       file_Dose1.open(outputDoseFileName1);
      file_Dose2.open(outPutDoseFileName2);
      file_Dose3.open(voxelPhantomInfo);
      
//       file_Dose<<"cellID"<<","<<"x"<<","<<"y"<<","<<"z"<<","<<"signalCountPrimaryMean"
//       <<","<<"signalCountPrimaryStd"<<endl; // file title
      
//       file_Dose1<<"cellID"<<","<<"x"<<","<<"y"<<","<<"z"<<","<<"signalCountPrimaryMean"
//       <<","<<"signalCountPrimaryStd"<<endl; // file title
      
      file_Dose2<<"cellID"<<","<<"x"<<","<<"y"<<","<<"z"<<","<<"signalCountPrimaryMean"
      <<","<<"signalCountPrimaryStd"<<endl; // file title
      
      file_Dose3<<"fNVoxelX"<<","<<"fNVoxelY"<<","<<"fNVoxelZ"<<","<<"fVoxelDimX"
      <<","<<"fVoxelDimY"<<","<<"fVoxelDimZ"<<endl;
      
	  
//       for (std::map<int, double>::iterator mitr_cell=cellSignalMeanMap.begin();mitr_cell!=cellSignalMeanMap.end();mitr_cell++)
//       {
// 
// 	G4double x = cellPositionXVec[mitr_cell->first] + doseTallVolumeISO_X;
// 	G4double y = cellPositionYVec[mitr_cell->first] + doseTallVolumeISO_Y;
// 	G4double z = cellPositionZVec[mitr_cell->first] + doseTallVolumeISO_Z;
// 	
// // 	cout<<"x = "<<x<<" y = "<<y<<"z = "<<z<<endl;
// 	
// 	file_Dose<<mitr_cell->first<<","<< x<<","<<y<<","<<z<<","<<cellSignalMeanMap[mitr_cell->first]<<","<<cellSignalStdMap[mitr_cell->first]<<endl;
//       }
      
      
      cout<<"size of cellPositionXVec1 "<<cellPositionXVec1.size()<<endl;
      cout<<"size of cellPositionXVec "<<cellPositionXVec.size()<<endl;
	 cout<<"size of cellSignalPrimaryMeanMap "<<cellSignalPrimaryMeanMap.size()<<endl;
	 cout<<"size of cellSignalMeanMap "<<cellSignalMeanMap.size()<<endl;
      
      
//       for (std::map<int, double>::iterator mitr_cell=cellSignalPrimaryMeanMap.begin();mitr_cell!=cellSignalPrimaryMeanMap.end();mitr_cell++)
//       {
// 
// 	G4double x = cellPositionXVec1[mitr_cell->first] + doseTallVolumeISO_X1;
// 	G4double y = cellPositionYVec1[mitr_cell->first] + doseTallVolumeISO_Y1;
// 	G4double z = cellPositionZVec1[mitr_cell->first] + doseTallVolumeISO_Z1;
// 	
// // 	cout<<"x = "<<x<<" y = "<<y<<"z = "<<z<<endl;
// 	
// 	file_Dose1<<mitr_cell->first<<","<< x<<","<<y<<","<<z<<","<<cellSignalPrimaryMeanMap[mitr_cell->first]<<","<<cellSignalPrimaryStdMap[mitr_cell->first]<<endl;
//       }
      
      
      
    G4int fNVoxelX = detector->GetVoxelNumberXYZ().getX();
    G4int fNVoxelY = detector->GetVoxelNumberXYZ().getY();
    G4int fNVoxelZ = detector->GetVoxelNumberXYZ().getZ();
    
    G4double fVoxelDimX = detector->GetVoxelDimension().getX();
    G4double fVoxelDimY = detector->GetVoxelDimension().getY();
    G4double fVoxelDimZ = detector->GetVoxelDimension().getZ();
    
    G4double fVoxelDimMinX = detector->GetVoxelDimMinXYZ().getX();
    G4double fVoxelDimMinY = detector->GetVoxelDimMinXYZ().getY();
    G4double fVoxelDimMinZ = detector->GetVoxelDimMinXYZ().getZ();
    
    G4double fVoxelDimMaxX = detector->GetVoxelDimMaxXYZ().getX();
    G4double fVoxelDimMaxY = detector->GetVoxelDimMaxXYZ().getY();
    G4double fVoxelDimMaxZ = detector->GetVoxelDimMaxXYZ().getZ();
    
    
    // add these lines to consider the phantom rotation to aline the phantom correctly in IEC 
    G4int fNVoxelY_tm = fNVoxelY; 
    G4int fNVoxelZ_tm = fNVoxelZ; 
    G4double fVoxelDimY_tm = fVoxelDimY;
    G4double fVoxelDimZ_tm = fVoxelDimZ;
    G4double fVoxelDimMinY_tm = fVoxelDimMinY;
    G4double fVoxelDimMinZ_tm = fVoxelDimMinZ;
    G4double fVoxelDimMaxY_tm = fVoxelDimMaxY;
    G4double fVoxelDimMaxZ_tm = fVoxelDimMaxZ;
    
    
    fNVoxelZ = fNVoxelY_tm;
    fNVoxelY = fNVoxelZ_tm;
    fVoxelDimY = fVoxelDimZ_tm;
    fVoxelDimZ = fVoxelDimY_tm;
    fVoxelDimMinY = fVoxelDimMinZ_tm;
    fVoxelDimMinZ = fVoxelDimMinY_tm;
    fVoxelDimMaxY = fVoxelDimMaxZ_tm;
    fVoxelDimMaxZ = fVoxelDimMaxY_tm;
    
    
    
    cout<<"fNVoxelX = "<<fNVoxelX <<" fNVoxelY = "<<fNVoxelY<<" fNVoxelZ = "<<fNVoxelZ<<endl;
    cout<<"fVoxelDimMinX = "<<fVoxelDimMinX<<" fVoxelDimMinY = "<<fVoxelDimMinY<<" fVoxelDimMinZ = "<<fVoxelDimMinZ<<endl;
    cout<<"fVoxelDimMaxX = "<<fVoxelDimMaxX<<" fVoxelDimMaxY = "<<fVoxelDimMaxY<<" fVoxelDimMaxZ = "<<fVoxelDimMaxZ<<endl;
    cout<<"fVoxelDimX = "<<fVoxelDimX<<" fVoxelDimY = "<<fVoxelDimY<<" fVoxelDimZ = "<<fVoxelDimZ<<endl;
//     exit(0);
    
    file_Dose3<<fNVoxelX<<","<<fNVoxelY<<","<<fNVoxelZ<<","<<fVoxelDimX
      <<","<<fVoxelDimY<<","<<fVoxelDimZ<<endl;
    
    for (std::map<G4int, G4double>::iterator mitr_cell = cellSignalScatterMeanMap.begin(); mitr_cell!=cellSignalScatterMeanMap.end(); mitr_cell++)
    {
      G4int copyNb = mitr_cell->first;
//       G4double x = fVoxelDimX*(copyNb%fNVoxelX);
//       G4double z = fVoxelDimY*((copyNb/fNVoxelX)%fNVoxelY);
//       G4double y = fVoxelDimZ*(copyNb/(fNVoxelX*fNVoxelY));
      
//         y[i] = fVoxelDimY*(i%fNVoxelY)
// 	x[i] = fVoxelDimX*((i/fNVoxelY)%fNVoxelX)
// 	z[i] = fVoxelDimZ*(i/(fNVoxelX*fNVoxelY))  
      
//       G4double y = fVoxelDimY*(copyNb%fNVoxelY); // uncomment these lines for printing x,y, z
//       G4double x = fVoxelDimX*((copyNb/fNVoxelY)%fNVoxelX);
//       G4double z = fVoxelDimZ*(copyNb/(fNVoxelX*fNVoxelY));
//       file_Dose2<<mitr_cell->first<<","<< x<<","<<y<<","<<z<<","<<cellSignalScatterMeanMap[mitr_cell->first]<<","<<cellSignalScatterStdMap[mitr_cell->first]<<endl; use this to print x, y, z
      file_Dose2<<mitr_cell->first<<","<<cellSignalScatterMeanMap[mitr_cell->first]<<","<<cellSignalScatterStdMap[mitr_cell->first]<<endl;
    }
    
//   z[i] = fVoxelDimZ*(i%fNVoxelZ)
//   x[i] = fVoxelDimX*((i/fNVoxelZ)%fNVoxelX)
//   y[i] = fVoxelDimY*(i/(fNVoxelX*fNVoxelZ))
      
      
      cout<<"Finished Analysis File Writing."<<endl;
      file_Dose.close();
      file_Dose1.close();
      file_Dose2.close();
      file_Dose3.close();
      
      
     myCollection.ClearCollectedSignalMap(); // delete the collected map after data processing
    }

}

void B3aRunAction::CollectEventCountMap(int eventID, std::map<int, double> sCurrentMap,std::map<int, double> cCurrentMap)
{
  eventSCurrentMap[eventID] = sCurrentMap;
  eventCCurrentMap[eventID] = cCurrentMap;
}

void B3aRunAction::SumDetectorEnergy(int detectorID, G4double energy, G4double energySquare)
{
//   fSumDetectorEnergy[detectorID] +=  energy;
}

pair< map< int, G4double >, map< int, G4double > > B3aRunAction::GetMCResults(map< int, G4double > countMap, map< int, G4double > countSquaredMap, int historyNumber)
{
  
     std::map<int, G4double > cellEnergyMap = countMap;
     std::map<int, G4double > cellEnergySquaredMap = countSquaredMap;
     int nofEvents = historyNumber;
  
      std::map<int,G4double>::iterator mitr_cell;
      std::map<int, double> cellSignalMeanMap;
      std::map<int, double> cellSignalStdMap;
     
     for (mitr_cell=cellEnergyMap.begin();mitr_cell!=cellEnergyMap.end();mitr_cell++)
     {
//        cout<<"cell energy is "<<mitr_cell->second<<endl;
//        cout<<"cell energy 2 is "<<cellEnergySquaredMap[mitr_cell->first]<<endl;
       cellSignalMeanMap[mitr_cell->first] = mitr_cell->second/nofEvents;
       G4double meanCountSquared = cellEnergySquaredMap[mitr_cell->first]/nofEvents;
       G4double std = sqrt(meanCountSquared-(cellSignalMeanMap[mitr_cell->first])*(cellSignalMeanMap[mitr_cell->first]));
       cellSignalStdMap[mitr_cell->first] = std/sqrt(nofEvents)*1.96;
//        cout<<"the mean energy is "<<cellSignalMeanMap[mitr_cell->first]<<endl;
//        cout<<"the mean energy squared is "<<meanCountSquared<<endl;
//        cout<<"std is "<<std<<endl;
//        exit(0);
     }
     
     std::pair<std::map<int, G4double>, std::map<int, G4double> > resultsMaps;
     resultsMaps.first = cellSignalMeanMap;
     resultsMaps.second = cellSignalStdMap;
     return resultsMaps;
     

}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
