
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// This the test function for the project

#include "RADCellSimulation.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include "G4RegionStore.hh"
#include <fstream>
#include <sstream>
#include <time.h>

#include "ArgumentInterpreter.hh"
#include "ReadLinacRTInformation.hh"
using namespace std;

// About usage of this function:
// we need to input two argeuments:
// The first one is simulation id, the second one is the control point filename

int main(int argc,char** argv)
{
//       for (int i = 0; i < argc; ++i) {
//         std::cout << argv[i] << std::endl;
//     }
    
    
    ArgumentInterpreter::SetIAEAPhaseSpaceFileIOType("read"); // could be read and write
    
    RADCellSimulation mySim;
    stringstream ss;
    ss<< argv[1];
    int i =0;
    ss>>i;
    cout<<"The phase space file index is: "<<i<<endl;
    stringstream ss1;
    ss1<<argv[2];
    G4String controlFileName;
    ss1>>controlFileName;
    cout<<controlFileName<<endl;
    
    cout<<"The controlFileName is "<<controlFileName<<endl;
//     exit(0);
//     cout<<controlFileName.substr(0,controlFileName.length()-5)<<endl;
    
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ControlPoints/";    
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ProsIMRT/";
    
    // this is for APPA plan
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/APPA/ControlPoints/";    
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/APPA/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/ProsIMRT/";
    
    
    // this is for prostate patient IMRT plan
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientIMRT/ControlPoints/";    
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientIMRT/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientIMRT/dicom/";
    
    
    // this is for brain patient IMRT plan
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/ControlPoints/";
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/dicom/";
    
    // this is for lung patinet IMRT plan
    
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientLung/ControlPoints/";
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientLung/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientLung/dicom/";
    
    
    // this is for prostate 002 IMRT plan
//        G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/";
    G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientProstate2/ControlPoints/";
    G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientProstate2/DoseOutput/";
    G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientProstate2/dicom/";
    
    // this is for fist abdomen patient
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientAbdomen1/ControlPoints/";
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientAbdomen1/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientAbdomen1/dicom/";
    
    
        // this is for brain patient IMRT plan
//     G4String pathOfControlPointsFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/ControlPoints/";
//     G4String pathOfDoseFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/DoseOutput/";
//     G4String pathofDicomFile = "/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/PatientBrainTumor2/dicom/";
    
    ArgumentInterpreter::SetConrolPointsFilePath(pathOfControlPointsFile);
    ArgumentInterpreter::SetDoseFilePath(pathOfDoseFile);
    ArgumentInterpreter::SetControlPointsFile(controlFileName);
    ArgumentInterpreter::SetDicomFilePath(pathofDicomFile);
    
    
    ReadLinacRTInformation myRead;
    
    myRead.ParseControlPointFile(controlFileName);
    
    cout<<"Mu: "<<myRead.GetMUInformation()<<endl;
    cout<<"sum Mu: "<<myRead.GetTotalMUInformation()<<endl;
    cout<<"GantryRtn: "<<myRead.GetGantryRtnInformation()<<endl;
    cout<<"CollRtn: "<<myRead.GetCollimationRtnInformation()<<endl;
    cout<<"CouchLat: "<<myRead.GetCoutchInformation().getX()<< " CouchVrt: "<<myRead.GetCoutchInformation().getZ()<<
    " CouchLng: "<<myRead.GetCoutchInformation().getY()<<endl;
    cout<<"ProxMLC: "<<endl;
    cout<<"The size of GetMLCDistInformation is "<<myRead.GetMLCDistInformation().size()<<endl;
    
    for (int i = 0;i< myRead.GetMLCProxInformation().size(); i++)
    {
      cout<<myRead.GetMLCProxInformation()[i]<<endl;
    }
    cout<<"DistMLC: "<<endl;
    for (int i = 0;i< myRead.GetMLCDistInformation().size(); i++)
    {
      cout<<myRead.GetMLCDistInformation()[i]<<endl;
    }
    
    
//        exit(0);
    
    ifstream ifile; // read  input file

    G4String radiationInputFileName;
    radiationInputFileName="run2.mac"; // input file 
    ifile.open(radiationInputFileName.c_str());
    
    stringstream ss_beamOn;
    stringstream ss2;
    string beamOnArgument;
    string inputFileName;
    
    std::string line;
    
    G4double totalMU = myRead.GetTotalMUInformation(); // get the total MU number
    G4double totalParticleNumber = 1E+9;// one billion particle default, you can change it 
    G4int particlePerMU = G4int(totalParticleNumber/totalMU);
    cout<<"particlePerMU = "<<particlePerMU<<endl;   
    ss_beamOn<<"/run/beamOn" <<" "<<G4int(particlePerMU*myRead.GetMUInformation()); // set up beam on partilce number
    
    ss2<<"run2_"<<"PNum"<<"_"<<G4int(particlePerMU*myRead.GetMUInformation())<<".mac";
    
    cout<<"simulation particle number for current control point is "<<G4int(particlePerMU*myRead.GetMUInformation())<<endl;
//     exit(0);
    
    inputFileName=ss2.str();
    ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file
    beamOnArgument=ss_beamOn.str();

    if (ifile.is_open())
    {
        while (!ifile.eof())
        {
            std::getline(ifile,line);
            std::istringstream ls(line);
            std::string firstString;
            std::string secondString;
            ls>>firstString>>secondString;
//                 cout<<"The first string is "<<firstString <<" The second string is "<<secondString<<endl;
            if (firstString=="/run/beamOn")
            {
                line=beamOnArgument;
            }
        
            ofile<<line<<endl;
        }
        
    }
    ifile.close();
    
//     exit(0);
    
    mySim.SetAndBuildDetectorArray(i); // the first source detector pair line simulation
    mySim.RADCellSimulationInitialize(argc, argv);
    mySim.EnergyDistributionCalculation("gui",inputFileName);
    
    
 
    
// //                 mySim.EnergyDistributionCalculation("gui",inputFileName,true); // run in gui mode
    
//     for (int i = 1; i<12; i++)
//     {
//       mySim.UpdateGeometryInitialize();
//       mySim.SetAndBuildDetectorArray(i*110);
//       mySim.UpdateGeometryFinalize();
//       mySim.EnergyDistributionCalculation("out",inputFileName);
//     }
// //     
//     mySim.EnergyDistributionCalculation("gui",inputFileName);
    
    return 0;
}

