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
  TChain* decimated = new TChain("headTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run < 257 || run > 263){
      TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%d.root", run, run);
      headChain->Add(fileName);

      fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);    
      decimated->Add(fileName);
    }
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);


  decimated->BuildIndex("eventNumber");


  const int numPulses=50;
  
  std::ofstream outFile("anita3OverwrittenEventInfo.txt", std::ofstream::out);
  outFile << "eventNumber\tfakeTreeEntry" << std::endl;
  
  Long64_t numMinBias = 0;
  // runing over events counts this many...
  Long64_t countedMinBiasEvents = 1813681;

  UInt_t seed = 123985; // mashed keyboard with hands
  // UInt_t seed = 0;
  TRandom3 rnd(seed);

  // don't print this number on the final go...
  const int N = rnd.Uniform(10, 15);
  // const int N = 10;
  
  // now divide minBiasEvents into N sections
  Long64_t numPerSection = countedMinBiasEvents / (N);

  std::vector<Int_t> indicesOfMinBiasEvents;
  std::vector<Int_t> fakeTreeEntries;
  std::vector<Int_t> fakeTreeEntriesAvailable(numPulses, 1);
  for(int i=0; i < N; i++){

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

  std::vector<Long64_t> boundaries;// {0, 10020905, 16529635, 25346226, 33950133, 41851692, 49780722, 56722096, 64172375, 71274663, 78630035, 78630542};
  
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  // ProgressBar p2(startEntry - maxEntry);
  for(Long64_t entry=startEntry; entry<maxEntry; entry++){
    headChain->GetEntry(entry);
    Int_t isMinBias = header->getTriggerBitSoftExt();
    if(isMinBias > 0){
      if(numMinBias % numPerSection == 0){
  	boundaries.push_back(entry);
  	// std::cerr << entry << std::endl;
      }
      numMinBias++;      
    }
    // p2.inc(startEntry, maxEntry);
  }

  boundaries.push_back(maxEntry);  
  // std::cout << " I got " << boundaries.size() << " boundaries" << std::endl << "{";
  // for(auto& e : boundaries){
  //   std::cout << e << ", ";
  // }
  // std::cout << "};" << std::endl;
  
  ProgressBar p(N);  
  for(Long64_t i=0; i < N; i++){

    bool selectedEventIsInDecimatedDataSet = true;
    while(selectedEventIsInDecimatedDataSet){
      numMinBias = 0;
      Long64_t indexOfMinBiasEvents = rnd.Uniform(0, numPerSection);
      
      for(Long64_t entry=boundaries.at(i); entry<boundaries.at(i+1); entry++){

	headChain->GetEntry(entry);

	Int_t isMinBias = header->getTriggerBitSoftExt();
	// std::cerr << entry << "\t" << boundaries.at(i) << "\t" << boundaries.at(i+1) << std::endl;
	if(isMinBias > 0){
	  if(indexOfMinBiasEvents==numMinBias){

	    Int_t isDec = decimated->GetEntryWithIndex(header->eventNumber);
	    if(isDec < 0){
	      outFile << header->eventNumber << "\t" << fakeTreeEntries.at(selectedEvents) << std::endl;
	      // std::cout << header->eventNumber << "\t" << header->run << "\t" << fakeTreeEntries.at(selectedEvents) << std::endl;
	      selectedEvents++;
	      selectedEventIsInDecimatedDataSet = false;
	      break;
	    }
	    else{
	      numMinBias = 0;
	      selectedEventIsInDecimatedDataSet = true;
	      break;
	    }
	  };
	  numMinBias++;
	  if(selectedEvents>=N){ // we're done here
	    break;
	  }
	}
      }
    }
    p.inc(i, N);
  }

  outFile.close();

  return 0;

}


