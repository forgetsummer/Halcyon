#include "ReadLinacRTInformation.hh"
#include <fstream>
#include <sstream>
#include "ArgumentInterpreter.hh"

void ReadLinacRTInformation::ParseControlPointFile(G4String fileName)
{
//   cout<<"get here"<<endl;
  fileName = ArgumentInterpreter::GetControlFilePath() + fileName;
  cout<<"fileName is "<<fileName<<endl;
  
  
  ifstream cpFile;
  cpFile.open(fileName);
  stringstream proximalDistance;
  stringstream distalDistance;
  string line;
  
  controlPointInfo.mu = 0;
  controlPointInfo.sum_mu = 0;
  controlPointInfo.gantryRtn = 0;
  controlPointInfo.collimationRtn = 0;
  controlPointInfo.coutchInfo = G4ThreeVector(0,0,0);
  controlPointInfo.proxMLCInfo.push_back(0);
  controlPointInfo.distMLCInfo.push_back(0);
  
  
  ControlPointsInformation cpInfo;
    if (cpFile.is_open())
    {
      while (getline(cpFile,line))
      {
// 	cout<<" cpInfo: "<<line<<endl;
	std::istringstream iss(line);
	std::vector<std::string> words;

	while (std::getline(iss, line, ' '))
	{
	  words.push_back(line);
	  
	}
// 	proximalDistance<<line;

// 	for (int i = 0; i<words.size(); i++)
// 	{
// 	  cout<<words[i]<<endl;
// 	}
	
	for (int i = 0; i<words.size(); i++)
	{
	  if (words[i] == "Mu")
	  {
	    stringstream smu;
	    stringstream ssum_mu;
	    smu<<words[i+1];
	    ssum_mu<<words[i+2];
	    G4double cmu;
	    G4double csum_mu;
	    smu>>cmu;
	    ssum_mu>>csum_mu;
	    cpInfo.mu =cmu;
	    cpInfo.sum_mu = csum_mu;
	  }
	  if (words[i] == "GantryRtn")
	  {
	    stringstream sgRtn;
	    sgRtn<<words[i+1];
	    G4double cgRtn;
	    sgRtn>>cgRtn;
	    cpInfo.gantryRtn = cgRtn;
	  }
	  if(words[i] == "CollRtn")
	  {
	    stringstream scRtn;
	    scRtn<<words[i+1];
	    G4double ccRtn;
	    scRtn>>ccRtn;
	    cpInfo.collimationRtn = ccRtn;
	  }
	  if (words[i] == "CouchLat")
	  {
	    stringstream schLat;
	    schLat<<words[i+1];
	    G4double x;
	    schLat>>x;
	    cpInfo.coutchInfo.setX(x);
	  }
	  if (words[i] == "CouchVrt")
	  {
	    stringstream schVrt;
	    schVrt<<words[i+1];
	    G4double x;
	    schVrt>>x;
	    cpInfo.coutchInfo.setZ(x);
	  }
	  if (words[i] == "CouchLng")
	  {
	    stringstream schLng;
	    schLng<<words[i+1];
	    G4double x;
	    schLng>>x;
	    cpInfo.coutchInfo.setY(x);
	  }
	  if (words[i] == "Prox")
	  {
	    stringstream ss( words[i+1] );
	    vector<string> result;
	    while( ss.good() )
	    {
		string substr;
		getline( ss, substr, ',' );
		stringstream lfpos;
		lfpos<<substr;
		G4double x;
		lfpos>>x;
		cpInfo.proxMLCInfo.push_back(x);
	    }
	  }
	  if (words[i] == "Dist")
	  {
	    stringstream ss( words[i+1] );
	    vector<string> result;
	    while( ss.good() )
	    {
		string substr;
		getline( ss, substr, ',' );
		stringstream lfpos;
		lfpos<<substr;
		G4double x;
		lfpos>>x;
		cpInfo.distMLCInfo.push_back(x);
	    }
	  }
	}

      }
    }
    cpFile.close();
    controlPointInfo = cpInfo;

}


G4double ReadLinacRTInformation::GetMUInformation()
{
  return controlPointInfo.mu;

}
G4double ReadLinacRTInformation::GetTotalMUInformation()
{
  return controlPointInfo.sum_mu;
}

G4double ReadLinacRTInformation::GetCollimationRtnInformation()
{
  return controlPointInfo.collimationRtn;

}
G4ThreeVector ReadLinacRTInformation::GetCoutchInformation()
{
  return controlPointInfo.coutchInfo;

}
G4double ReadLinacRTInformation::GetGantryRtnInformation()
{
  return controlPointInfo.gantryRtn;
}

std::vector<G4double> ReadLinacRTInformation::GetMLCDistInformation()
{
  return controlPointInfo.distMLCInfo;
}
std::vector< G4double > ReadLinacRTInformation::GetMLCProxInformation()
{
  return controlPointInfo.proxMLCInfo;
}



