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
// $Id: B3StackingAction.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B3StackingAction.hh
/// \brief Definition of the B3StackingAction class

#ifndef B3StackingAction_h
#define B3StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4Region.hh"

/// Stacking action class : manage the newly generated particles
///
/// One wishes do not track secondary neutrino.Therefore one kills it 
/// immediately, before created particles will  put in a stack.

class B3StackingAction : public G4UserStackingAction
{
  public:
    B3StackingAction();
    virtual ~B3StackingAction();
     
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);   
  private:
  G4Region* fPhantomRegion; 
  
//       G4int               fKillSecondary;


    G4int               fPhotoGamma;
    G4int               fComptGamma;
//     G4int               fPhotoAuger;
//     G4int               fComptAuger;
//     G4int               fPixeGamma;
//     G4int               fPixeAuger;

//     G4int               fElectronDNAGamma;
//     G4int               fElectronDNAAuger;
//     G4int               fProtonDNAGamma;
//     G4int               fProtonDNAAuger;
//     G4int               fHydrogenDNAGamma;
//     G4int               fHydrogenDNAAuger;
//     G4int               fAlphaDNAGamma;
//     G4int               fAlphaDNAAuger;
//     G4int               fAlphaPlusDNAGamma;
//     G4int               fAlphaPlusDNAAuger;
//     G4int               fHeliumDNAGamma;
//     G4int               fHeliumDNAAuger;
//     G4int               fGenericIonDNAGamma;
//     G4int               fGenericIonDNAAuger;

    G4bool              fIDdefined;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

