// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015
// This is a class for running the radiation transportation in cells

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "B3DetectorConstruction.hh"
#include "LinacSimulationDetectorConstruction.hh"
#include "volumeConstruction.hh" // this one is for building the linac head according to the example: G4Linac-head
#include "MRLinacHead.hh" // this is for building the linac head according to the example: medical-linac in Geant4
#include "VarianClinac2100Head.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4Types.hh"
#include <string>

using namespace std;


class RADCellSimulation 
 {
      public: 
      RADCellSimulation();// the constructor
      ~RADCellSimulation(); // the destructor
      void SetAndBuildDetectorArray(int sourceDetectorPairID);
      void RADCellSimulationInitialize(int argc, char** argv); // function initializing the radiation transport simulation
      void RADCellSimulationInitializePyWrapper(int argc,const std::vector<std::string> & _vec_str); // for python wrapping RADCellSimulationInitialize()
      void EnergyDistributionCalculation(string runArgument,string inputFile); // function for running the radiation transport, the arguments are the running mode.
      void UpdateGeometryInitialize(); // the function initializing the process of changing the cell geometry between different runs
      void UpdateGeometryFinalize(); // the function finalizing the process of changing the cell geometry between different runs

      
      
      private:
      LinacSimulationDetectorConstruction* linacDetector;
      volumeConstruction* detector;
      MRLinacHead* detector1;
      VarianClinac2100Head* detector2;
      G4VisManager* visManager;
        #ifdef G4MULTITHREADED
        G4MTRunManager* runManager;
        #else
        G4RunManager* runManager;
        #endif
      G4UImanager* UImanager;
      G4UIExecutive* ui; 

      
      

 };
	
