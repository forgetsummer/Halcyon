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


#ifndef VarianClinac2100Head_h
#define VarianClinac2100Head_h 1

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


class VarianClinac2100Head
{
public:
	VarianClinac2100Head();
	~VarianClinac2100Head();
	
	void GetLinacHeadLogicVolume(G4VPhysicalVolume* physicalVolume);
	G4ThreeVector GetTargetPosition() {return targetPosition;};

	void reset();

	inline void setIsoCentre(G4double val){isoCentre=val;}
	inline void setidEnergy(G4int val){idEnergy=val;}

	inline int getidEnergy(){return idEnergy;}
	G4double getBeforeJaws_Z_PhaseSpacePosition(){return 215.;}
	void writeInfo();
private:
	G4double jaw1XAperture, jaw2XAperture, jaw1YAperture, jaw2YAperture, isoCentre; 
	inline void setJaw1X(G4double val){jaw1XAperture=val;}
	inline void setJaw2X(G4double val){jaw2XAperture=val;}
	inline void setJaw1Y(G4double val){jaw1YAperture=val;}
	inline void setJaw2Y(G4double val){jaw2YAperture=val;}
	
	inline void setLeavesAx(G4double val){leavesA.push_back(val);}
	inline void setLeavesBx(G4double val){leavesB.push_back(val);}
	
	std::vector <G4double> leavesA, leavesB;
	G4int idEnergy;

	G4Material * otherMaterials(const G4String materialName);
	void SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4double aperture, G4RotationMatrix *cRotation);
	void target();
	void primaryCollimator();
	void BeWindow();
	void flatteningFilter();
	void ionizationChamber();
	void mirror();
	void Jaw1X();
	void Jaw2X();
	void Jaw1Y();
	void Jaw2Y();
	void MLC();
	
	void PhaseSpacePlane(G4double zPosition); // add a virtual plane for showing the position of the phase space plane
	
	G4ThreeVector relativePositionOfLinacHead;

	G4VPhysicalVolume * PVWorld;
        G4VPhysicalVolume *targetA_phys;
        G4VPhysicalVolume *targetB_phys;
        G4VPhysicalVolume *UpperCollimator_phys;
        G4VPhysicalVolume *CylMinusCone_phys;
        G4VPhysicalVolume *BeWTubePV;
        G4VPhysicalVolume *FFL1A_1PV;
        G4VPhysicalVolume *FFL2_1PV;
        G4VPhysicalVolume *PCUtubeW1PV;
        G4VPhysicalVolume *PCUtubeP1PV;
        G4VPhysicalVolume *PCUtubeW2PV;
        G4VPhysicalVolume *PCUtubeP2PV;
        G4VPhysicalVolume *PCUtubeW3PV;
        G4VPhysicalVolume *PCUtubeP3PV;
        G4VPhysicalVolume *MirrorTubePV;
        G4VPhysicalVolume *phVol1X;
        G4VPhysicalVolume *phVol2X;
        G4VPhysicalVolume *phVol1Y;
        G4VPhysicalVolume *phVol2Y;
        G4VPhysicalVolume *leafPhys;
	
	G4ThreeVector targetPosition;

};

#endif

