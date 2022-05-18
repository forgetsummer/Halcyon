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
// $Id: B3StackingAction.cc 66536 2012-12-19 14:32:36Z ihrivnac $
// 
/// \file B3StackingAction.cc
/// \brief Implementation of the B3StackingAction class

#include "B3StackingAction.hh"
#include "G4Track.hh"
#include "G4NeutrinoE.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4VSolid.hh"
#include "G4RegionStore.hh"
#include <iostream>
#include <stdlib.h>
#include "G4PhysicsModelCatalog.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3StackingAction::B3StackingAction():fPhantomRegion(0),fPhotoGamma(-1),fComptGamma(-1)

{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3StackingAction::~B3StackingAction()
{ 
  fPhantomRegion = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
B3StackingAction::ClassifyNewTrack(const G4Track* track)
{

  //keep primary particle
  if (track->GetParentID() == 0) return fUrgent;
  
//   if(!fIDdefined) {
//     fIDdefined = true;
// //     fPhotoGamma = G4PhysicsModelCatalog::GetIndex("phot_fluo");
// //     fComptGamma = G4PhysicsModelCatalog::GetIndex("compt_fluo");
// //     fPhotoAuger = G4PhysicsModelCatalog::GetIndex("phot_auger");
// //     fComptAuger = G4PhysicsModelCatalog::GetIndex("compt_auger");
// //     fPixeGamma = G4PhysicsModelCatalog::GetIndex("gammaPIXE");
// //     fPixeAuger = G4PhysicsModelCatalog::GetIndex("e-PIXE");
// 
// 
//   }
  
//     G4int entries = G4PhysicsModelCatalog::Entries();
//     for (G4int i = 0; i < entries; ++i) {
//       cout << i << ": " << G4PhysicsModelCatalog::GetModelName(i) << endl;
//     }
//   G4int idx = track->GetCreatorModelID();

//   cout<<"fPhotoGamma is "<<fPhotoGamma<<endl;
//   cout<<"fComptGamma is "<<fComptGamma<<endl;
   

//     cout<< "idx is "<<idx<< " The creat model name: " << G4PhysicsModelCatalog::GetModelName(24)<<endl;

  //kill secondary neutrino
//   if (track->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
//   else return fUrgent;
  
  // kill the track if there is a scattering happen when the primary particle is inside the patient volume
  fPhantomRegion = G4RegionStore::GetInstance()->GetRegion("Phantom"); // get the phantom fPantomRegion
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  if(particleDefinition == G4Electron::Definition() ||
    particleDefinition == G4Gamma::Definition())
  {

      if(fPhantomRegion == 0)
	// target region is initialized after
	//detector construction instantiation
      {
	  G4Exception("TrackingAction","fPhantomRegion == 0",
		      FatalException,"No phantom has been found");
      }
      
      const G4ThreeVector& position = track->GetPosition(); // get the track position
      int N =  fPhantomRegion->GetNumberOfRootVolumes();
      std::vector<G4LogicalVolume*>::iterator it_logicalVolumeInRegion =
	    fPhantomRegion->GetRootLogicalVolumeIterator();
      bool inside_target = false;
      
      for(int i = 0; i < N ; i++, it_logicalVolumeInRegion++)
      {
	  EInside test_status = (*it_logicalVolumeInRegion)->GetSolid()->
	      Inside(position) ;
	  if(test_status == kInside)
	  {
	      inside_target = true;
	      break;
	  }
	  /*
	  else if (test_status == kSurface)
	  {
	  }
	  */
      }
      if (inside_target = true) // if the particle is inside phantom
      {
// 	cout<<"the particle is inside phantom "<<endl;
// 	exit(0);
// 	if (track->GetParentID()>0) // it is a secondary particle
// 	{
// // 	    cout<<"Kill the secondary particles inside phantom "<<endl;
// // 	    exit(0);
// // 	  return fKill;
// 	}
// 	if (track->GetCreatorModelName()=="Rayl")
// 	if (track->GetCreatorModelName()=="compt" || track->GetCreatorModelName() == "Rayl")
// 	{
// 	  return fKill;
// 	}
// 	else
// 	  return fUrgent;

      }
      
  }
  
//   cout<< "idx is "<<track->GetCreatorModelName()<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
