// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1
#include "G4Region.hh"
#include "G4UserSteppingAction.hh"

class B3aEventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SteppingAction : public G4UserSteppingAction
{
public:

  SteppingAction(B3aEventAction* eventAction);
  ~SteppingAction();
  
  virtual void UserSteppingAction(const G4Step*);
private:
  B3aEventAction* fEventAction;

  
};
#endif
