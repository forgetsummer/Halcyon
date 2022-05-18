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
// $Id: B3PrimaryGeneratorAction.cc 73744 2013-09-09 20:25:07Z asaim $
//
/// \file B3PrimaryGeneratorAction.cc
/// \brief Implementation of the B3PrimaryGeneratorAction class

#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh" // for using constants, such as pi here

#include "G4IAEAphspReader.hh"

#include "G4AutoLock.hh" // this is for thread safety of using multithreading mode for reading phase file by G4IAEAphspReader

#include "ArgumentInterpreter.hh"

#include "TruebeamHead.hh"

#include "ReadLinacRTInformation.hh"


using namespace CLHEP;
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


namespace { G4Mutex myLowEPrimGenMutex	= G4MUTEX_INITIALIZER;}
G4IAEAphspReader* B3PrimaryGeneratorAction::theIAEAReader = 0;

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
                    = particleTable->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(1*eV);    
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  
//   G4String fileName = "MRLinacPhaseSpaceFile_750_1BParticle";
  
  stringstream ss;
//   ss<<"/home/administrator/Geant4Projects/TruebeamSimulation/MRIdian-build/VarianTrueBeamPS6MVFFF/";
  ss<<"/home/administrator/Geant4Projects/TruebeamSimulation/HalcyonSimulation/Halcyon-build/VarianHalcyoonPS6MVFFF/";
  
  std::vector<G4String>fileNameVec;
  fileNameVec.push_back("TrueBeam_v2_6FFF_00");
  fileNameVec.push_back("TrueBeam_v2_6FFF_01");
  fileNameVec.push_back("TrueBeam_v2_6FFF_02");
  fileNameVec.push_back("TrueBeam_v2_6FFF_03");
  fileNameVec.push_back("TrueBeam_v2_6FFF_04");
  fileNameVec.push_back("TrueBeam_v2_6FFF_05");
  fileNameVec.push_back("TrueBeam_v2_6FFF_06");
  fileNameVec.push_back("TrueBeam_v2_6FFF_07");
  fileNameVec.push_back("TrueBeam_v2_6FFF_08");
  fileNameVec.push_back("TrueBeam_v2_6FFF_09");
  fileNameVec.push_back("TrueBeam_v2_6FFF_10");
  fileNameVec.push_back("TrueBeam_v2_6FFF_11");
  fileNameVec.push_back("TrueBeam_v2_6FFF_12");
  fileNameVec.push_back("TrueBeam_v2_6FFF_13");
  fileNameVec.push_back("TrueBeam_v2_6FFF_14");
  fileNameVec.push_back("TrueBeam_v2_6FFF_15");
  fileNameVec.push_back("TrueBeam_v2_6FFF_16");
  fileNameVec.push_back("TrueBeam_v2_6FFF_17");
  fileNameVec.push_back("TrueBeam_v2_6FFF_18");
  fileNameVec.push_back("TrueBeam_v2_6FFF_19");
  fileNameVec.push_back("TrueBeam_v2_6FFF_20");
  fileNameVec.push_back("TrueBeam_v2_6FFF_21");
  fileNameVec.push_back("TrueBeam_v2_6FFF_22");
  fileNameVec.push_back("TrueBeam_v2_6FFF_23");
  fileNameVec.push_back("TrueBeam_v2_6FFF_24");
  fileNameVec.push_back("TrueBeam_v2_6FFF_25");
  fileNameVec.push_back("TrueBeam_v2_6FFF_26");
  fileNameVec.push_back("TrueBeam_v2_6FFF_27");
  fileNameVec.push_back("TrueBeam_v2_6FFF_28");
  fileNameVec.push_back("TrueBeam_v2_6FFF_29");
  fileNameVec.push_back("TrueBeam_v2_6FFF_30");
  fileNameVec.push_back("TrueBeam_v2_6FFF_31");
  fileNameVec.push_back("TrueBeam_v2_6FFF_32");
  fileNameVec.push_back("TrueBeam_v2_6FFF_33");
  fileNameVec.push_back("TrueBeam_v2_6FFF_34");
  fileNameVec.push_back("TrueBeam_v2_6FFF_35");
  fileNameVec.push_back("TrueBeam_v2_6FFF_36");
  fileNameVec.push_back("TrueBeam_v2_6FFF_37");
  fileNameVec.push_back("TrueBeam_v2_6FFF_38");
  fileNameVec.push_back("TrueBeam_v2_6FFF_39");
  fileNameVec.push_back("TrueBeam_v2_6FFF_40");
  fileNameVec.push_back("TrueBeam_v2_6FFF_41");
  fileNameVec.push_back("TrueBeam_v2_6FFF_42");
  fileNameVec.push_back("TrueBeam_v2_6FFF_43");
  fileNameVec.push_back("TrueBeam_v2_6FFF_44");
  fileNameVec.push_back("TrueBeam_v2_6FFF_45");
  fileNameVec.push_back("TrueBeam_v2_6FFF_46");
  fileNameVec.push_back("TrueBeam_v2_6FFF_47");
  fileNameVec.push_back("TrueBeam_v2_6FFF_48");
  fileNameVec.push_back("TrueBeam_v2_6FFF_49");
  fileNameVec.push_back("TrueBeam_v2_6FFF_50");
  fileNameVec.push_back("TrueBeam_v2_6FFF_51");
  fileNameVec.push_back("TrueBeam_v2_6FFF_52");
  fileNameVec.push_back("TrueBeam_v2_6FFF_53");
  fileNameVec.push_back("TrueBeam_v2_6FFF_54");
  
  
  std::vector<G4String>fileNameVecHalcyon;
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_1");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_2");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_3");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_4");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_5");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_6");
  fileNameVecHalcyon.push_back("Halcyon_FFF_6MV_7");
  
  int nn = ArgumentInterpreter::GetSourceDetectorPairID();
  
  G4cout<<"nn = "<<nn<<G4endl;
//   ss<<"TrueBeam_v2_6X_01";
//   ss<<fileNameVec[nn];
//   ss<<"OF_1e7_phasespace_field";
//   ss<<"test_10k_field";
//   ss<<"1e8part_nophantom_field";
//   ss<<"Halycon_FFF_6MV_1";
  ss<<fileNameVecHalcyon[nn];
  
  G4String fileName = ss.str();
  
//   G4String fileName = "test_10k_BS1_field";
  cout<<"fileName is "<<fileName<<endl;
//   exit(0);
//   G4String fileName = "MRLinacPhaseSpaceFile_440_0.3";
//   G4String fileName = "TrueBeam6MV_02_28June2011_0";
  
  if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "read") // if we read the IAEA phace space file
  {
      G4AutoLock lock(&myLowEPrimGenMutex);
      if(!theIAEAReader) theIAEAReader = new G4IAEAphspReader(fileName);
      
//       theIAEAReader->SetGlobalPhspTranslation(G4ThreeVector(0,0,300)); // this makes that the phase space plane is 267mm from source
      
    G4String fileNameLinac =  ArgumentInterpreter::GetControlPointsFile(); // get linac information file name
    ReadLinacRTInformation readHlcy;
    readHlcy.ParseControlPointFile(fileNameLinac);// parse the file information

    G4double gantryRtn = readHlcy.GetGantryRtnInformation()*deg - 180*deg;
    G4double collRtn = readHlcy.GetCollimationRtnInformation()*deg ;
    cout<<"GantryRtn: "<<gantryRtn<<endl;
    cout<<"CollRtn: " << collRtn<<endl;
    
    G4ThreeVector oldPosition = G4ThreeVector(0,0, 733*mm); // original position of the phase space plane
    G4ThreeVector oldPositionReflect = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn, collRtn );
//     cout<<"xr = "<<oldPositionReflect.getX()<<" yr = "<<oldPositionReflect.getY()<<" zr = "<<oldPositionReflect.getZ()<<endl;
    
    G4ThreeVector newPosition = GetCenterPositionOfRotatedMLC(oldPosition,gantryRtn + pi, collRtn + pi );
//     cout<<"xr = "<<newPosition.getX()<<" yr = "<<newPosition.getY()<<" zr = "<<newPosition.getZ()<<endl;
//     exit(0);
    G4double x0 = newPosition.getX()- oldPositionReflect.getX();
    G4double y0 = newPosition.getY()- oldPositionReflect.getY();
    G4double z0 = newPosition.getZ()- oldPositionReflect.getZ();
    
//     cout<<"x0 = "<<x0<<" y0= "<<y0<<" z0= "<<z0<<endl;
    
//     exit(0);
    G4ThreeVector globalXYZ = G4ThreeVector(x0,y0,z0);
      theIAEAReader->SetGlobalPhspTranslation(globalXYZ); // 
      theIAEAReader->SetRotationOrder(321);
//       theIAEAReader->SetRotationY(180*deg);
      
      theIAEAReader->SetRotationY(gantryRtn);
      theIAEAReader->SetRotationZ(-1*collRtn);
//       
      
//       theIAEAReader->SetTotalParallelRuns(16); // 16 fragments
//       theIAEAReader->SetParallelRun(7); // 8rd fragment of the PSF, for debugging 
/*      
      theIAEAReader->SetTotalParallelRuns(16); // 8 fragments
      theIAEAReader->SetParallelRun(7); // 3rd fragment of the PSF*/
      theIAEAReader->SetTimesRecycled(24);
      
  }
  
//       G4AutoLock lock(&myLowEPrimGenMutex);
//       if(!theIAEAReader) theIAEAReader = new G4IAEAphspReader(fileName);
  
  
  
  
//   G4ThreeVector psfShift(0., 0.,-50*cm);//lespace de phase est 
//   theIAEAReader->SetGlobalPhspTranslation(psfShift);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun;
   G4AutoLock lock(&myLowEPrimGenMutex);
  if(theIAEAReader) 
  {
    delete theIAEAReader; theIAEAReader = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
     G4int n_particle = 1;
//    fParticleGun  = new G4ParticleGun(n_particle);
 
   // default particle kinematic
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4String particleName;
//    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-"); // specifying the particle type of electron gun
   G4ParticleDefinition* particle = particleTable->FindParticle(particleName="chargedgeantino");  
   fParticleGun->SetParticleDefinition(particle);
   
   
   
   
   if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "write")
   {
      G4double mu_z = -1; //
      G4double mu_x = 0;
      G4double mu_y = 0;
      std::vector<std::vector<double> > transitionMatrix;
      transitionMatrix =   MatrixInitialization(3,3);
      G4double theta = atan(10/2/90.0);
      G4String jawName = "Jaw2X";
      if (jawName == "Jaw1X")
	{
	  
	  transitionMatrix[0][0] = cos(theta);
	  transitionMatrix[0][1] = 0;
	  transitionMatrix[0][2] = sin(theta);
	  transitionMatrix[1][0] = 0;
	  transitionMatrix[1][1] = 1;
	  transitionMatrix[1][2] = 0;
	  transitionMatrix[2][0] = -sin(theta);
	  transitionMatrix[2][1] = 0;
	  transitionMatrix[2][2] = cos(theta);

	}
	if (jawName == "Jaw2X")
	{
	  transitionMatrix[0][0] = cos(theta);
	  transitionMatrix[0][1] = 0;
	  transitionMatrix[0][2] = -sin(theta);
	  transitionMatrix[1][0] = 0;
	  transitionMatrix[1][1] = 1;
	  transitionMatrix[1][2] = 0;
	  transitionMatrix[2][0] = sin(theta);
	  transitionMatrix[2][1] = 0;
	  transitionMatrix[2][2] = cos(theta);
	  
	}
	if (jawName == "Jaw1Y")
	{
	  transitionMatrix[0][0] = 1;
	  transitionMatrix[0][1] = 0;
	  transitionMatrix[0][2] = 0;
	  transitionMatrix[1][0] = 0;
	  transitionMatrix[1][1] = cos(theta);
	  transitionMatrix[1][2] = sin(theta);
	  transitionMatrix[2][0] = 0;
	  transitionMatrix[2][1] = -sin(theta);
	  transitionMatrix[2][2] = cos(theta);
	  
	  
	}
	if (jawName == "Jaw2Y")
	{
	  transitionMatrix[0][0] = 1;
	  transitionMatrix[0][1] = 0;
	  transitionMatrix[0][2] = 0;
	  transitionMatrix[1][0] = 0;
	  transitionMatrix[1][1] = cos(theta);
	  transitionMatrix[1][2] = -sin(theta);
	  transitionMatrix[2][0] = 0;
	  transitionMatrix[2][1] = sin(theta);
	  transitionMatrix[2][2] = cos(theta);
	  
	}
      
      
      std::vector<std::vector<double> > oldPosition;
      oldPosition = MatrixInitialization(3,1);
      G4ThreeVector initialPosition = G4ThreeVector(mu_x,mu_y,mu_z);
      
      oldPosition[0][0] = initialPosition.getX();
      oldPosition[1][0] = initialPosition.getY();
      oldPosition[2][0] = initialPosition.getZ();
      
      std::vector<std::vector<double> > newPosition;
      newPosition = MatrixMultiplication(transitionMatrix,oldPosition);
      
      G4ThreeVector positionAfterRotation;
      positionAfterRotation.setX(newPosition[0][0]);
      positionAfterRotation.setY(newPosition[1][0]);
      positionAfterRotation.setZ(newPosition[2][0]);
//       cout<<"ui = "<<positionAfterRotation.getX()<<endl;
//       cout<<"uj = "<<positionAfterRotation.getY()<<endl;
//       cout<<"uk = "<<positionAfterRotation.getZ()<<endl;
//       exit(0);
      
      if (jawName == "Jaw1X" || jawName == "Jaw2X")
      {
	positionAfterRotation.setX(-newPosition[0][0] + 2*newPosition[0][0]*G4UniformRand());
      }
      
      if (jawName == "Jaw1Y" || jawName == "Jaw2Y")
      {
	positionAfterRotation.setY(-newPosition[1][0] + 2*newPosition[1][0]*G4UniformRand());
      }
//       fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mu_x,mu_y,mu_z)); // set the photon shooting angle
      fParticleGun->SetParticleMomentumDirection(positionAfterRotation);

	G4double meanKineticEnergy = 6*MeV;
	G4double stdDev = 0.3*MeV;
	G4double kineticEnergy = RandGauss::shoot(meanKineticEnergy, stdDev);
    //     cout<<"kineticEnergy is "<<kineticEnergy<<endl;
	fParticleGun->SetParticleEnergy(kineticEnergy); // set the photon energy, here we use a X-ray spectrum 

      G4double X_source;
      G4double Y_source;
      G4double Z_source = 90*cm;
      G4double planeSourceRadius = 0*cm;
      G4double phi = 2*pi*G4UniformRand();
      G4double radius = planeSourceRadius*G4UniformRand();
      X_source = radius*cos(phi);
      Y_source = radius*sin(phi);
      
      fParticleGun->SetParticlePosition(G4ThreeVector(X_source,Y_source,Z_source)); // set the position of the point source at the center
		
      //create vertex
      //
      fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  
    if (ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType() == "read")
    {
        // using phase space file to generate particle
      
         
      if(theIAEAReader)
      {
	G4AutoLock lock(&myLowEPrimGenMutex);
	theIAEAReader->GeneratePrimaryVertex(anEvent);
      }
    }
    

}

std::vector< std::vector< double > > B3PrimaryGeneratorAction::MatrixInitialization(int m, int n)
{
  
    std::vector<std::vector<double> > M;
    for (int i=0;i<m;i++ )
    {
        M.push_back(std::vector<double>());
        for (int j=0;j<n;j++)
        {
            M[i].push_back(0); // initial result
        }
    }
    return M;
}

std::vector< std::vector< double > > B3PrimaryGeneratorAction::MatrixMultiplication(std::vector< std::vector< double > > M1, std::vector< std::vector< double > > M2)
{
    std::vector<std::vector<double> >  result;
    int m = M1.size();
    int n = M2[0].size();
    result = MatrixInitialization(m,n);
//     cout<<"m = "<<m<<" n = "<<n<<endl;
//     PrintMatrix(result);
    
    for (int i = 0;i<M1.size();i++)
    {
        for (int j = 0;j< M2[0].size();j++)
        {
            for (int k=0;k<M1[0].size();k++)
            {
                result[i][j] += M1[i][k] * M2[k][j];
            }
        }
    }
    return result;  
}

G4ThreeVector B3PrimaryGeneratorAction::GetCenterPositionOfRotatedMLC(G4ThreeVector initialPosition, G4double gantryRtn, G4double collRtn)
{
  std::vector<std::vector<double> > transitionMatrixRg; // transition matrix of gantry angle
  std::vector<std::vector<double> > transitionMatrixRc; // transition matrix of collimation angle
  transitionMatrixRg =   MatrixInitialization(3,3);
  transitionMatrixRc =   MatrixInitialization(3,3);
    
//     transitionMatrixRg[0][0] = 1;
//     transitionMatrixRg[0][1] = 0;
//     transitionMatrixRg[0][2] = 0;
//     transitionMatrixRg[1][0] = 0;
//     transitionMatrixRg[1][1] = cos(gantryRtn);
//     transitionMatrixRg[1][2] = -sin(gantryRtn);
//     transitionMatrixRg[2][0] = 0;
//     transitionMatrixRg[2][1] = sin(gantryRtn);
//     transitionMatrixRg[2][2] = cos(gantryRtn); // gantry rotates with X
    
        
    transitionMatrixRg[0][0] = cos(gantryRtn);
    transitionMatrixRg[0][1] = 0;
    transitionMatrixRg[0][2] = sin(gantryRtn);
    transitionMatrixRg[1][0] = 0;
    transitionMatrixRg[1][1] = 1;
    transitionMatrixRg[1][2] = 0;
    transitionMatrixRg[2][0] = -sin(gantryRtn);
    transitionMatrixRg[2][1] = 0;
    transitionMatrixRg[2][2] = cos(gantryRtn); // gantry rotates with Y
    
    
    
    transitionMatrixRc[0][0] = cos(collRtn);
    transitionMatrixRc[0][1] = -sin(collRtn);
    transitionMatrixRc[0][2] = 0;
    transitionMatrixRc[1][0] = sin(collRtn);
    transitionMatrixRc[1][1] = cos(collRtn);
    transitionMatrixRc[1][2] = 0;
    transitionMatrixRc[2][0] = 0;
    transitionMatrixRc[2][1] = 0;
    transitionMatrixRc[2][2] = 1; // collimation rotates with Z
    
  
//   PrintMatrix(transitionMatrix);
  
    std::vector<std::vector<double> > transitionMatrixRg0; // transition matrix of gantry angle
    std::vector<std::vector<double> > transitionMatrixRc0; // transition matrix of collimation angle
    transitionMatrixRg0 =   MatrixInitialization(3,3);
    transitionMatrixRc0 =   MatrixInitialization(3,3);
    G4double gantryRtn0 = -pi;
//     transitionMatrixRg0[0][0] = 1;
//     transitionMatrixRg0[0][1] = 0;
//     transitionMatrixRg0[0][2] = 0;
//     transitionMatrixRg0[1][0] = 0;
//     transitionMatrixRg0[1][1] = cos(gantryRtn0);
//     transitionMatrixRg0[1][2] = -sin(gantryRtn0);
//     transitionMatrixRg0[2][0] = 0;
//     transitionMatrixRg0[2][1] = sin(gantryRtn0);
//     transitionMatrixRg0[2][2] = cos(gantryRtn0);
    
    
    
    transitionMatrixRg0[0][0] = cos(gantryRtn0);
    transitionMatrixRg0[0][1] = 0;
    transitionMatrixRg0[0][2] = sin(gantryRtn0);
    transitionMatrixRg0[1][0] = 0;
    transitionMatrixRg0[1][1] = 1;
    transitionMatrixRg0[1][2] = 0;
    transitionMatrixRg0[2][0] = -sin(gantryRtn0);
    transitionMatrixRg0[2][1] = 0;
    transitionMatrixRg0[2][2] = cos(gantryRtn0); // gantry rotates with Y
    
    G4double collRtn0 = -pi;
    
    transitionMatrixRc0[0][0] = cos(collRtn0);
    transitionMatrixRc0[0][1] = -sin(collRtn0);
    transitionMatrixRc0[0][2] = 0;
    transitionMatrixRc0[1][0] = sin(collRtn0);
    transitionMatrixRc0[1][1] = cos(collRtn0);
    transitionMatrixRc0[1][2] = 0;
    transitionMatrixRc0[2][0] = 0;
    transitionMatrixRc0[2][1] = 0;
    transitionMatrixRc0[2][2] = 1;
    
  
  std::vector<std::vector<double> > oldPosition;
  oldPosition = MatrixInitialization(3,1);
  oldPosition[0][0] = initialPosition.getX();
  oldPosition[1][0] = initialPosition.getY();
  oldPosition[2][0] = initialPosition.getZ();
  
//   PrintMatrix(oldPosition);
  
  
  
//   oldPosition = MatrixMultiplication(transitionMatrixRc0,oldPosition); // restore to original position
//   oldPosition = MatrixMultiplication(transitionMatrixRg0,oldPosition); // restore to original position
  
  std::vector<std::vector<double> > newPosition;
  newPosition = MatrixMultiplication(transitionMatrixRc,oldPosition); //implement collimation rotation
  newPosition = MatrixMultiplication(transitionMatrixRg,newPosition); // implement gantry rotation rotation
  
  
  G4ThreeVector positionAfterRotation;
  positionAfterRotation.setX(newPosition[0][0]);
  positionAfterRotation.setY(newPosition[1][0]);
  positionAfterRotation.setZ(newPosition[2][0]);
  
  return positionAfterRotation;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

