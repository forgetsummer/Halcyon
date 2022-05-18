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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef TruebeamHead_h
#define TruebeamHead_h 1

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ProductionCuts.hh"
#include "G4GenericTrap.hh"

class TruebeamHead
{
public:
	TruebeamHead();
	~TruebeamHead();
	
	void GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume);
	G4ThreeVector GetTargetPosition() {return targetPosition;};



private:
 
  
	void DefineMaterials();
	G4double jaw1XAperture, jaw2XAperture, jaw1YAperture, jaw2YAperture, isoCentre; 
	inline void setJaw1X(G4double val){jaw1XAperture=val;}
	inline void setJaw2X(G4double val){jaw2XAperture=val;}
	inline void setJaw1Y(G4double val){jaw1YAperture=val;}
	inline void setJaw2Y(G4double val){jaw2YAperture=val;}
	
	void MLC();
	G4double GetMLCPositionA(G4double kx,G4double ssd, G4double R,  G4double y);
	G4double GetMLCPositionB(G4double kx,G4double ssd, G4double R,  G4double y);

	G4ThreeVector GetCenterPositionOfRotatedJaw(G4String jawName, G4ThreeVector initialPosition, G4double theta); // theta is rotation angle
	std::vector<std::vector<double> > MatrixInitialization(int m, int n);
	std::vector<std::vector<double> > MatrixMultiplication(std::vector<std::vector<double> > M1, std::vector<std::vector<double> > M2); // doing matrix multiplication
	void PrintMatrix(std::vector<std::vector<double> > M);
	
	void Jaw1X();
	void Jaw2X();
	void Jaw1Y();
	void Jaw2Y();
	
	void SetFieldSize(G4double fieldSize);
	void SetSSD(G4double ssd);
	G4double fieldSize;
	G4double SSD;
	
	void BasePlate();
	
	G4ThreeVector GetCenterPositionOfRotatedMLC(G4ThreeVector initialPosition, G4double gantryRtn, G4double collRtn);

	
	void PhaseSpacePlane(G4double zPosition); // add a virtual plane for showing the position of the phase space plane
	
	void SourceMLCLeakBlock(G4double zPosition);
	
	G4double MLCPosCorrection(G4double x);
	
	G4ThreeVector relativePositionOfLinacHead;

	G4VPhysicalVolume * PVWorld;

        G4VPhysicalVolume *phVol1X;
        G4VPhysicalVolume *phVol2X;
        G4VPhysicalVolume *phVol1Y;
        G4VPhysicalVolume *phVol2Y;
	
	G4ThreeVector targetPosition;

};

#endif

