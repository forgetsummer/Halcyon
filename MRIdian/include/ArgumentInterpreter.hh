// 
// This is code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015

#ifndef ArgumentInterpreter_h
#define ArgumentInterpreter_h 1
#include "G4String.hh"

struct stat;class ArgumentInterpreter
{
    public:
        static void SetArgument(G4String arg);
	static void SetSourceDetectorPairID(int pairID);
        static G4String GetOutPutFileName();
	static G4int GetSourceDetectorPairID();
	
	static void SetIAEAPhaseSpaceFileIOType(G4String readOrWrite);
	static G4String GetIAEAPhaseSpaceFileIOType();
	
	static void SetControlPointsFile(G4String controlFile);
	static G4String GetControlPointsFile();
	
	static void SetConrolPointsFilePath(G4String fpathControlFile);
	static G4String GetControlFilePath();
	
	static void SetDoseFilePath(G4String fpathDoseFile);
	static G4String GetDoseFilePath();
	
	static void SetDicomFilePath(G4String fpathDicomFile);
	static G4String GetDicomFilePath();
private:
        static G4String argument;
	static G4int sourceDetectorPairID;
	static G4String writeOrReadIAEAPhaseSpaceFile;
	static G4String ControlPointsFile;
	static G4String pathControlFile;
	static G4String pathDoseFile;
	static G4String pathDicomFile;

        
};
#endif