// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to select WAIS pulses, swap their polarizations info and put them into a new tree
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "FancyFFTs.h"

int main(int argc, char* argv[]){

  // Runs near WAIS divide 
  // const Int_t firstRun = 331;
  // const Int_t lastRun = 354;
  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);  

  //*************************************************************************
  // Set up input
  //*************************************************************************  

  TChain* headChain = new TChain("headTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.unblind.root", run, run);
    headChain->Add(fileName);

  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  const int numPulses=50;


  
  std::ofstream outFile("anita3OverwrittenEventInfo.txt", std::ofstream::out);
  outFile << "eventNumber\tfakeTreeEntry" << std::endl;
  
  Long64_t numMinBias = 0;
  // runing over events counts this many...
  Long64_t countedMinBiasEvents = 1813681;

  UInt_t seed = 123984; // mashed keyboard with hands
  TRandom3 rnd(seed);

  // don't print this number on the final go...
  const int N = rnd.Uniform(10, 15);
  // const int N = 10;
  
  // now divide minBiasEvents into N sections
  const int numPerSection = countedMinBiasEvents / (N);

  std::vector<Int_t> indicesOfMinBiasEvents;
  std::vector<Int_t> fakeTreeEntries;
  std::vector<Int_t> fakeTreeEntriesAvailable(numPulses, 1);
  for(int i=0; i < N; i++){

    Int_t j = rnd.Uniform(0, numPerSection);
    Int_t minBiasIndex = numPerSection*i + j;
    indicesOfMinBiasEvents.push_back(minBiasIndex);

    Int_t fakeTreeEntry = -1;
    bool unusedEntry=false;
    while(unusedEntry==false){
      Int_t tryThisEntry = rnd.Uniform(0, numPulses);
      if(fakeTreeEntriesAvailable.at(tryThisEntry) == 1){
	fakeTreeEntry = tryThisEntry;
	fakeTreeEntriesAvailable.at(tryThisEntry) = 0;
	unusedEntry = true;
      }
    }
    fakeTreeEntries.push_back(fakeTreeEntry);
  }
  Int_t selectedEvents = 0;


  // This variable tracks which entry the replaced event will map to
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry=startEntry; entry<maxEntry; entry++){

    headChain->GetEntry(entry);

    Int_t isMinBias = header->getTriggerBitSoftExt();

    if(isMinBias > 0){
      if(indicesOfMinBiasEvents.at(selectedEvents)==numMinBias){
	// std::cerr << numMinBias << "\t" << header->eventNumber << "\t" << header->run << std::endl;
	// std::cerr << numMinBias << "\t" << header->eventNumber << "\t"
	// 	  << fakeTreeEntries.at(selectedEvents) << "\t"
	// 	  << header->run << std::endl;
	
	outFile << header->eventNumber << "\t" << fakeTreeEntries.at(selectedEvents) << std::endl;
	selectedEvents++;
      };
      numMinBias++;
      if(selectedEvents>=N){ // we're done here
	entry=maxEntry-1;
      }
    }
	
    p.inc(entry, maxEntry);
  }

  outFile.close();
  std::cout << "I counted " << numMinBias << " minimum bias triggered events " << std::endl;

  return 0;

}


