
#ifndef ReadLinacRTInformation_h
#define ReadLinacRTInformation_h 1
#include <vector>
#include <G4String.hh>
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

using namespace std;


struct ControlPointsInformation
{
  G4double mu;// this stores the MU number for current control point
  G4double sum_mu;//this stores the total MU number in the plan
  std::vector<G4double> proxMLCInfo;
  std::vector<G4double> distMLCInfo;
  G4double gantryRtn;
  G4double collimationRtn;
  G4ThreeVector coutchInfo;
};

class ReadLinacRTInformation
{
public: 
    void ParseControlPointFile(G4String fileName); 
    G4double GetMUInformation();
    G4double GetTotalMUInformation();
    std::vector<G4double> GetMLCProxInformation();
    std::vector<G4double> GetMLCDistInformation();
    G4double GetGantryRtnInformation();
    G4double GetCollimationRtnInformation();
    G4ThreeVector GetCoutchInformation();
private:
  ControlPointsInformation controlPointInfo;
  
    
  
};




#endif