#ifndef projectDataTypeGlobals_h
#define projectDataTypeGlobals_h 1

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <vector>
#include <map>
#include <utility> 
using namespace std;

struct EDEPInformation
{
  G4double x;
  G4double y;
  G4double z;
  G4double edep;
};



#endif
