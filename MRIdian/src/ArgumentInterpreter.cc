#include "ArgumentInterpreter.hh"
#include <sstream>
#include <vector>
using namespace std;

G4String ArgumentInterpreter::argument;// Initialize the argument here
G4int ArgumentInterpreter::sourceDetectorPairID;
G4String ArgumentInterpreter::writeOrReadIAEAPhaseSpaceFile;
G4String ArgumentInterpreter::ControlPointsFile;
G4String ArgumentInterpreter::pathControlFile;
G4String ArgumentInterpreter::pathDoseFile;
G4String ArgumentInterpreter::pathDicomFile;


void ArgumentInterpreter::SetArgument(G4String arg)
{
    argument = arg;
}

void ArgumentInterpreter::SetSourceDetectorPairID(int pairID)
{
  sourceDetectorPairID = pairID;
}

G4String ArgumentInterpreter::GetOutPutFileName()
{
    G4String argument_case;
    std::vector<G4String>argument_vector;
    stringstream ss(argument);
    while (ss>>argument_case)
    {
        argument_vector.push_back(argument_case);
    }
    
    G4String filename = argument_vector[0];

    return filename;

}
G4int ArgumentInterpreter::GetSourceDetectorPairID()
{
  return sourceDetectorPairID;
}

void ArgumentInterpreter::SetIAEAPhaseSpaceFileIOType(G4String readOrWrite)
{
  writeOrReadIAEAPhaseSpaceFile = readOrWrite;

}
G4String ArgumentInterpreter::GetIAEAPhaseSpaceFileIOType()
{
  return writeOrReadIAEAPhaseSpaceFile;

}

void ArgumentInterpreter::SetControlPointsFile(G4String controlFile)
{
  ControlPointsFile = controlFile;

}

G4String ArgumentInterpreter::GetControlPointsFile()
{
  return ControlPointsFile;
}

void ArgumentInterpreter::SetConrolPointsFilePath(G4String fpathControlFile)
{
  pathControlFile = fpathControlFile;
}
G4String ArgumentInterpreter::GetControlFilePath()
{
  return pathControlFile;
}

void ArgumentInterpreter::SetDoseFilePath(G4String fpathDoseFile)
{
  pathDoseFile = fpathDoseFile;
}
G4String ArgumentInterpreter::GetDoseFilePath()
{
  return pathDoseFile;
}

void ArgumentInterpreter::SetDicomFilePath(G4String fpathDicomFile)
{
  pathDicomFile = fpathDicomFile;
}
G4String ArgumentInterpreter::GetDicomFilePath()
{
  return pathDicomFile;
}





