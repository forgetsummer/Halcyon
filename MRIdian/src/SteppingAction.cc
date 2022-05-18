
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class



#include "SteppingAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include <vector>
#include <sstream>
#include "G4EventManager.hh" // using this class to get the current event ID
#include "G4Track.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "B3aEventAction.hh"
#include "projectDataTypeGlobals.hh"

#include "G4IAEAphspWriter.hh"

#include "ArgumentInterpreter.hh"
#include "G4TransportationManager.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(B3aEventAction* eventAction)
  :G4UserSteppingAction(),fEventAction(eventAction)/// We should use construct initializer
{
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
//  delete fPhantomRegion;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  
//   G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetPushVerbosity(0);// set up the verbosity level of G4Navigator class
  
    G4double flagParticle=0.;
  G4double flagProcess=0.;
  G4double x,y,z,xp,yp,zp;


  const G4String& processName = step->GetPostStepPoint()->
      GetProcessDefinedStep()->GetProcessName();
 /////////////////////////////////////////////////////////////////////////////////
 /// testing the distribution of phase space file

   if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "write")
   {
       G4IAEAphspWriter::GetInstance()->UserSteppingAction(step);
   } 
  
  
  

  
//       // get volume of the current step
//       G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//       G4String edepVolume = volume->GetName();
//       
// //       cout<<"the volume is "<<edepVolume<< endl;
      if (step->GetTrack()->GetGlobalTime()/second >10)
      {
	  cout<<"the tracking time is "<<step->GetTrack()->GetGlobalTime()<<endl;
      }
      
//       if (edepVolume == "phantom")
//       {
// 	  const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
// 	  G4int ix = touchable->GetReplicaNumber(1);
// 	  G4int iy = touchable->GetReplicaNumber(2);
// 	  G4int iz = touchable->GetReplicaNumber(0);
// 	  cout<<"ix = "<<ix<<" iy = "<<iy<<" iz = "<<iz<<endl;
// 	  
// 	    G4StepPoint* thePostPoint = step->GetPostStepPoint();
// 	    G4StepPoint* thePrePoint = step->GetPreStepPoint();
// 	    x=step->GetPreStepPoint()->GetPosition().x()/mm;
// 	    y=step->GetPreStepPoint()->GetPosition().y()/mm;
// 	    z=step->GetPreStepPoint()->GetPosition().z()/mm;
// 	    xp=step->GetPostStepPoint()->GetPosition().x()/mm;
// 	    yp=step->GetPostStepPoint()->GetPosition().y()/mm;
// 	    zp=step->GetPostStepPoint()->GetPosition().z()/mm;
// 	    
// 	    cout<<"x= "<<x<<" y = "<<y<<" z= "<<z<<endl;
// 	    cout<<"xp= "<<xp<<" yp= "<<yp<<" zp=s "<<zp<<endl;
// // 	  G4int tmp = fNy;
// // 	  if (tmp) return iy*fNx*fNz+ix*fNz+iz;
// // 	  else return iy*fNx*fNz+ix*fNz+iz;
// 	
//       }
      
      

// // 	  if (edepVolume=="phantom" )
//          if (edepVolume=="phantom" || edepVolume == "unionV") // the physical volume name of dicom voxel is phantom
// 	 {
// 	    G4StepPoint* thePostPoint = step->GetPostStepPoint();
// 	    G4StepPoint* thePrePoint = step->GetPreStepPoint();
// 	    x=step->GetPreStepPoint()->GetPosition().x()/mm;
// 	    y=step->GetPreStepPoint()->GetPosition().y()/mm;
// 	    z=step->GetPreStepPoint()->GetPosition().z()/mm;
// 	    xp=step->GetPostStepPoint()->GetPosition().x()/mm;
// 	    yp=step->GetPostStepPoint()->GetPosition().y()/mm;
// 	    zp=step->GetPostStepPoint()->GetPosition().z()/mm;



// 	    
// 	    
// // 	    cout<<"the track material is "<<step->GetTrack()->GetMaterial()->GetName()<<endl;
// 	    
// 	
// // 	cout<<"x= "<<x<<" y = "<<y<<" z= "<<z<<endl;
// // 	cout<<"xp= "<<xp<<" yp= "<<yp<<" zp=s "<<zp<<endl;
// 	
// 	    G4double stepLength = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp));
// 	    if (stepLength<1e-5)
// 	    {
// 	      step->GetTrack()->SetTrackStatus(fStopAndKill); // when the step is not moving, we kill it
// 	    }
// 	    if (step->GetTrack()->GetCurrentStepNumber()>10000) //here we use this to monitor the step number in a track
// 	    {
// 	      cout<<"the current step number in this track is "<<step->GetTrack()->GetCurrentStepNumber()<<endl;
// 	    }
// 	    
// 	 }
  
  
   /////////////////////////////////////////////////////////////////////////////////////
   ////testing the phase space file spatial distribution
//    if( particleDefinition == G4Gamma::Definition())
//    {
// 
//     // get volume of the current step
//       G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//       G4String edepVolume = volume->GetName();
//       
//       if (edepVolume=="waterTank_phys")
//       {
// 	
// 	x=step->GetPreStepPoint()->GetPosition().x()/mm;
// 	y=step->GetPreStepPoint()->GetPosition().y()/mm;
// 	z=step->GetPreStepPoint()->GetPosition().z()/mm;
// 	xp=step->GetPostStepPoint()->GetPosition().x()/mm;
// 	yp=step->GetPostStepPoint()->GetPosition().y()/mm;
// 	zp=step->GetPostStepPoint()->GetPosition().z()/mm;
// 	
// 	  // energy deposit
// 	G4double edep = step->GetTotalEnergyDeposit()/eV;// get deposited energy in eV
// 	EDEPInformation edepEventInfo;
// 	edepEventInfo.x=xp;
// 	edepEventInfo.y=yp;
// 	edepEventInfo.z=zp;
// 	edepEventInfo.edep=edep;
// 	fEventAction->CollectEDEPInformation(edepEventInfo);
// 	
// 	
//       }
//     
//     }
    
  //////////////////////////////////////////////////////////////////////////////  
//  if (fEventAction->GetCurrentEventID()>7000)
//  {
//    cout<<fEventAction->GetCurrentEventID()<<endl;
//    
//    	x=step->GetPreStepPoint()->GetPosition().x()/mm;
// 	y=step->GetPreStepPoint()->GetPosition().y()/mm;
// 	z=step->GetPreStepPoint()->GetPosition().z()/mm;
// 	xp=step->GetPostStepPoint()->GetPosition().x()/mm;
// 	yp=step->GetPostStepPoint()->GetPosition().y()/mm;
// 	zp=step->GetPostStepPoint()->GetPosition().z()/mm;
// 	
// 	cout<<"x= "<<x<<" y= "<<y<<" z= "<<z<<endl;
// 	cout<<"xp= "<<xp<<" yp= "<<yp<<" zp= "<<zp<<endl;
// 	cout<<"step energy "<<step->GetPostStepPoint()->GetTotalEnergy()<<endl;
// 	cout<<"deposited energy "<<step->GetTotalEnergyDeposit()/eV<<endl;
// //    cout<<processName<<endl;
// //    if (processName.empty())
// //    {
// //      cout<<processName<<endl;
// // //      cout<<""<<endl;
// // //      cerr<<"what is going on "<<endl;
// //      
// //   }
// //    cerr<<"what is going on "<<endl;
// // //    cout<<processName<<endl;
//  }

//********************below is used to calculate the dose distribution for the water phantom
   
//   if (processName!="Transportation") // supposed to decide by checking processName, but there is some unknown bug for reading the halcyron phase space file
//   if (step->GetTotalEnergyDeposit()>0) // we decide by checkinig the deposited energy, if bigger than zero, we do the followings
//   {
// 
// 
// 
//   // get volume of the current step
//     G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//   
//   // energy deposit
//     G4double edep = step->GetTotalEnergyDeposit()/eV;// get deposited energy in eV
//     
//     G4String edepVolume = volume->GetName();
// 
//     
//     if (edepVolume=="waterTank_phys")
//     {
//         
//       x=step->GetPreStepPoint()->GetPosition().x()/mm;
//       y=step->GetPreStepPoint()->GetPosition().y()/mm;
//       z=step->GetPreStepPoint()->GetPosition().z()/mm;
//       xp=step->GetPostStepPoint()->GetPosition().x()/mm;
//       yp=step->GetPostStepPoint()->GetPosition().y()/mm;
//       zp=step->GetPostStepPoint()->GetPosition().z()/mm;
// //       cout<<"the volume is "<<edepVolume<<endl;
// //       cout<<"x = "<<x<<" "<<"y = "<<y<<" z = "<<z <<endl;
//       
// //       cout<<"The step length is "<<step->GetStepLength()/mm<< " mm"<<endl;
//       EDEPInformation edepEventInfo;
//       edepEventInfo.x=xp;
//       edepEventInfo.y=yp;
//       edepEventInfo.z=zp;
//       edepEventInfo.edep=edep;
//       fEventAction->CollectEDEPInformation(edepEventInfo);
//       
//     }
//     
//   }

}
