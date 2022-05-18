// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 

#include "RADCellSimulation.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "B3aActionInitialization.hh"
#include "B3PhysicsList.hh"
#include "G4UIQt.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4Material.hh"
#include <sstream>
#include <iostream>
#include "ArgumentInterpreter.hh"
#include "G4SDManager.hh"
using namespace std;

RADCellSimulation::RADCellSimulation()
{
    linacDetector = new LinacSimulationDetectorConstruction();
    detector = new volumeConstruction();
    detector1 = new MRLinacHead();
    detector2 = new VarianClinac2100Head();
}


RADCellSimulation::~RADCellSimulation()
{
    delete visManager;
    delete runManager;
//     delete ui;

}


void RADCellSimulation::SetAndBuildDetectorArray(int sourceDetectorPairID)
{
//   detector->BuildAndSetDetectorArray(sourceDetectorPairID);
  string filename;
  stringstream ss;
  ss<<sourceDetectorPairID;
  filename = "sourceDetectorPairID_"+ss.str();
  ArgumentInterpreter::SetArgument(filename); // set the file name for output in each run
  ArgumentInterpreter::SetSourceDetectorPairID(sourceDetectorPairID); // set the sourceDetectorPairID in each run
  
  cout<<"FINISH SetAndBuildDetectorArray"<<endl;
}

void RADCellSimulation::RADCellSimulationInitializePyWrapper(int argc,const std::vector<std::string> & _vec_str)
{
    cerr<<"argc="<<argc<<endl;
    cerr<<"_vec_str.size()="<<_vec_str.size()<<endl;
    char** argv  = (char**)&_vec_str[0];
    cerr<<"RADCellSimulationInitializePyWrapper argv  [0]="<<argv[0]<<endl;
    
    for (int i = 0 ; i < _vec_str.size() ; ++i){
        cerr<<"argument["<<i<<"]="<<argv[i]<<endl;
    }
    
    RADCellSimulationInitialize(argc,argv);  
}
//----------------------------------------------------------------------
// Below is for function of EnergyDistributionCalculation()
//----------------------------------------------------------------------



void RADCellSimulation::RADCellSimulationInitialize(int argc, char** argv)
{
    // Choose the Random engine
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
//   //choose the Random engine
//   CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
//   //set random seed with system time
//   G4long seed = time(NULL);
//   CLHEP::HepRandom::setTheSeed(seed);

  #ifdef G4MULTITHREADED
   // G4MTRunManager* runManager = new G4MTRunManager;
    runManager= new G4MTRunManager;
    runManager->SetNumberOfThreads(40);// set up the number of number of threads for multithread simulation
  #else
   // G4RunManager* runManager = new G4RunManager;
    runManager= new G4RunManager;
  #endif
  
//   runManager->SetUserInitialization(detector2);//here we use the internal variable detector in class 
  runManager->SetUserInitialization(linacDetector);//here we use the internal variable detector in class 


  runManager->SetUserInitialization(new B3PhysicsList);

  // User action initialization
  
  runManager->SetUserInitialization(new B3aActionInitialization);
  


  // Initialize G4 kernel
  
    runManager->Initialize();
 
  // Initialize visualization 
  visManager = new G4VisExecutive;

  visManager->Initialize();  // comment this if you want to run radcell in CC3D
  ui = new G4UIExecutive(argc, argv); // 
  UImanager = G4UImanager::GetUIpointer();  // Get the pointer to the User Interface manager
//   UImanager->ApplyCommand("/control/execute vis.mac");

}

void RADCellSimulation::UpdateGeometryInitialize()
{
   runManager->ReinitializeGeometry(true);
//    G4GeometryManager::GetInstance()->OpenGeometry();
//    G4RegionStore::GetInstance()->Clean();// delete the regions in store
  
   cout<<"FINISHED UpdateGeometryInitialize"<<endl;
   
}
void RADCellSimulation::UpdateGeometryFinalize()
{
   
  runManager->DefineWorldVolume(detector->Construct());
  runManager->PhysicsHasBeenModified();
  runManager->GeometryHasBeenModified();
//   delete G4SDManager::GetSDMpointer();
//   G4SDManager::GetSDMpointer()->Activate("crystal",false);
//   G4SDManager::GetSDMpointer()->Activate("patient",false);
  cout<<"FINISHED UpdateGeometryFinalize"<<endl;

}


void RADCellSimulation::EnergyDistributionCalculation(string runArgument,string inputFile)
{
// Process macro or start UI session
  
    stringstream argument(runArgument); // get run argument
    G4String argument_case;
    std::vector<G4String>argument_vector;
    stringstream commandArgument;
    commandArgument<<"/control/execute" <<" "<<inputFile;
    string runCommand;
    runCommand=commandArgument.str();
    
  
    
    while (argument>>argument_case)
    {
        argument_vector.push_back(argument_case);
       
    }

    

  //
    if (argument_vector[0]=="gui")// running in gui interactive mode 
    {
          // Define UI terminal for interactive mode  
        //G4UIExecutive* ui = new G4UIExecutive(argc, argv);

	UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart(); 
       delete ui;
    }
    if (argument_vector[0]=="out")// running in output file mode
    { 
      cout<<"Before ApplyCommand"<<endl;
        UImanager->ApplyCommand(runCommand);
        
//         UImanager->ApplyCommand("/control/execute microdosimetry.in");
        
    }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
