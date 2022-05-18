#ifndef G4RunCollectTriggeredSignals_h
#define G4RunCollectTriggeredSignals_h 1
#include <map>
#include "G4UnitsTable.hh"
#include "G4Types.hh"
#include "G4Cache.hh"

using namespace std;


class G4RunCollectTriggeredSignals
{
public:
  static void CollectTriggeredSignal(int cellID, G4double energy);
  static void CollectTriggeredSignalPrimary(int cellID, G4double energy);
  static void CollectTriggeredSignalScatter(int cellID, G4double energy);
  std::map<int, G4double > GetWholeRunSignalCountMap () {return wholeRunSignalCount;}
  std::map<int, G4double > GetWholeRunSignalCountSquaredMap (){return wholeRunSignalCountSquared;}
  std::map<int, G4double > GetWholeRunSignalCountPrimaryMap(){return wholeRunSignalCountPrimary;}
  std::map<int, G4double > GetWholeRunSignalCountPrimarySquaredMap(){return wholeRunSignalCountPrimarySquared;}
  std::map<int, G4double > GetWholeRunSignalCountScatterMap(){return wholeRunSignalCountScatter;}
  std::map<int, G4double > GetWholeRunSignalCountScatterSquaredMap(){return wholeRunSignalCountScatterSquared;}
  

  
  
  void ClearCollectedSignalMap();
private:
  
  static  std::map<int, G4double > wholeRunSignalCount;
  static  std::map<int, G4double > wholeRunSignalCountSquared;
  static  std::map<int, G4double > wholeRunSignalCountPrimary;
  static  std::map<int, G4double > wholeRunSignalCountPrimarySquared;
  static  std::map<int, G4double > wholeRunSignalCountScatter;
  static  std::map<int, G4double > wholeRunSignalCountScatterSquared;
  
  
};

#endif