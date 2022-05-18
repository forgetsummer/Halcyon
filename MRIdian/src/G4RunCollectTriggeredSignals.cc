#include "G4RunCollectTriggeredSignals.hh"
namespace { G4Mutex myLowEPrimGenMutex	= G4MUTEX_INITIALIZER;}


// EventDoseTallyMap G4RunDoseAnalysis::doseTallyMap; // declare an instance of this static variable
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCount;
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCountSquared;
 
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCountPrimary;
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCountPrimarySquared;
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCountScatter;
  std::map<int, G4double > G4RunCollectTriggeredSignals::wholeRunSignalCountScatterSquared;


 
 void G4RunCollectTriggeredSignals::CollectTriggeredSignal(int cellID, G4double energy)
{
  G4AutoLock lock(&myLowEPrimGenMutex); // here we try to add a lock mechanism to do the non-thread safe operation. by RL 9/24/2019
  wholeRunSignalCount[cellID] = wholeRunSignalCount[cellID] + energy;
  lock.lock();
  wholeRunSignalCountSquared[cellID] = wholeRunSignalCountSquared[cellID] + energy*energy;
}

void G4RunCollectTriggeredSignals::ClearCollectedSignalMap()
{
  wholeRunSignalCount.clear();
  wholeRunSignalCountSquared.clear();
  wholeRunSignalCountPrimary.clear();
  wholeRunSignalCountPrimarySquared.clear();
  wholeRunSignalCountScatter.clear();
  wholeRunSignalCountScatterSquared.clear();
};

void G4RunCollectTriggeredSignals::CollectTriggeredSignalPrimary(int cellID, G4double energy)
{
  G4AutoLock lock(&myLowEPrimGenMutex);
  wholeRunSignalCountPrimary[cellID] =   wholeRunSignalCountPrimary[cellID] + energy;
  lock.lock();
  wholeRunSignalCountPrimarySquared[cellID] =   wholeRunSignalCountPrimarySquared[cellID] + energy*energy;

}
void G4RunCollectTriggeredSignals::CollectTriggeredSignalScatter(int cellID, G4double energy)
{
  G4AutoLock lock(&myLowEPrimGenMutex);
  wholeRunSignalCountScatter[cellID] = wholeRunSignalCountScatter[cellID]+ energy;
  lock.lock();
  wholeRunSignalCountScatterSquared[cellID] =  wholeRunSignalCountScatterSquared[cellID] + energy*energy;

}

