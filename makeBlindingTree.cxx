// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to find peak cross correlation offsets between antenna pairs in pulses from Wais Divide.
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"

int main()
{

  // Runs near WAIS divide 
  // const Int_t firstRun = 331;
  // const Int_t lastRun = 354;
  const Int_t firstRun = 352;
  const Int_t lastRun = 352;  

  
  TChain* calEventChain = new TChain("eventTree");
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);
  calEventChain->BuildIndex("eventNumber");
  
  TString outFileName = "blindEventTreeFile.root";
  TFile* outFile = new TFile(outFileName, "recreate");
  TTree* blindEventTree = new TTree("blindEventTree", "blindEventTree");
  UInt_t hashOfEventNumberToBeReplaced = 0;
  CalibratedAnitaEvent* changedCalEvent = NULL;
  blindEventTree->Branch("hashOfEventNumberToBeReplaced", &hashOfEventNumberToBeReplaced);
  blindEventTree->Branch("event", &changedCalEvent);

  
  std::ifstream blindFile("listOfEventsToSwap.txt");
  char firstLine[1024];
  blindFile.getline(firstLine,1024);
  std::cout << "the top line reads: " << firstLine << std::endl;

  std::vector<UInt_t> eventNumberToBeReplaced;
  std::vector<UInt_t> eventNumberToReplaceItWith;
  UInt_t replaceMe = 0;
  UInt_t withMe = 0;
  while(blindFile >> replaceMe >> withMe) {
    eventNumberToBeReplaced.push_back(replaceMe);
    eventNumberToReplaceItWith.push_back(withMe);
  }

  //  calEventChain->Show(0);
  
  ProgressBar p(eventNumberToBeReplaced.size());
  for(UInt_t entry = 0; entry < eventNumberToBeReplaced.size(); entry++){
    replaceMe = eventNumberToBeReplaced.at(entry);
    withMe = eventNumberToReplaceItWith.at(entry);

    TString replaceMeString = TString::Format("%u", replaceMe);
    hashOfEventNumberToBeReplaced = replaceMeString.Hash();
    calEventChain->GetEntryWithIndex(withMe);

    changedCalEvent = (CalibratedAnitaEvent*) calEvent->Clone();

    for(Int_t polInd=0; polInd<AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t inputPol = (AnitaPol::AnitaPol_t) polInd;
      // AnitaPol::AnitaPol_t outputPol = inputPol == AnitaPol::kHorizontal ? AnitaPol::kVertical : AnitaPol::kHorizontal;
      AnitaPol::AnitaPol_t outputPol = inputPol;
      
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	Int_t inputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, inputPol);
	Int_t outputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, outputPol);

	// std::cout << inputPol << "\t" << inputIndex << "\t\t"
	// 	  << outputPol << "\t" << outputIndex << "\t\t"
	// 	  << inputIndex/9 << "\t" << outputIndex/9 << std::endl;

	for(int samp=0; samp<260; samp++){
	  changedCalEvent->data[outputIndex][samp] = calEvent->data[inputIndex][samp];
	}
	changedCalEvent->xMax[outputIndex] = calEvent->xMax[inputIndex];
	changedCalEvent->xMin[outputIndex] = calEvent->xMin[inputIndex];
	changedCalEvent->mean[outputIndex] = calEvent->mean[inputIndex];
	changedCalEvent->rms[outputIndex] = calEvent->rms[inputIndex];
      }
    }
    
    blindEventTree->Fill();
  
    p++;
  }
  outFile->Write();
  //  blindEventTree->Show(0);
  
  outFile->Close();

  return 0;
}
