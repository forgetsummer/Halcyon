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
// $Id: B3aRunAction.hh 99559 2016-09-27 07:02:21Z gcosmo $
//
/// \file B3aRunAction.hh
/// \brief Definition of the B3aRunAction class

#ifndef B3aRunAction_h
#define B3aRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <map>
#include <utility>      // std::pair, std::make_pair
#include "CellLayoutInitializer.hh"

using namespace std;

/// Run action class

class B3aRunAction : public G4UserRunAction
{
  public:
    B3aRunAction();
    virtual ~B3aRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void CountEvent()           { fGoodEvents += 1; };
    void SumDose(G4double dose) { fSumDose += dose; };  
    
    void SumDetectorEnergy(int detectorID, G4double energy, G4double energySquare);
    
    void CollectEventCountMap(int eventID, std::map<int, double> sCurrentMap,std::map<int, double> cCurrentMap);
    


private:
    G4Accumulable<G4int>    fGoodEvents;
    G4Accumulable<G4double> fSumDose;  

    std::map<int, G4Accumulable<G4double>  > fSumDetectorEnergy;
    std::map<int, G4Accumulable<G4double>  > fSumDetectorEnergySquare;
    
    
    std::map<int, std::map<int, double> > eventSCurrentMap; //for storing the flat surface current
    //here we define a two dimensional map for storing the current of each detector in each event. The first key is event ID
    // the second key is the cell ID
    
    std::map<int, std::map<int, double> > eventCCurrentMap; // for storing the cell current 
    
    std::pair<std::map<int,G4double>, std::map<int, G4double> > GetMCResults(std::map<int,G4double> countMap, std::map<int,G4double> countSquaredMap, int historyNumber);
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

