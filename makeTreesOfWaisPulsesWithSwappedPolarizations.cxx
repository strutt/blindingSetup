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

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"

int main(int argc, char* argv[]){

  // Runs near WAIS divide 
  // const Int_t firstRun = 331;
  // const Int_t lastRun = 354;
  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " [run]" << std::endl;
  }
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun;

  //*************************************************************************
  // Set up input
  //*************************************************************************  

  TChain* calEventChain = new TChain("eventTree");
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
    
  }
  CalibratedAnitaEvent* calEventIn = NULL;
  calEventChain->SetBranchAddress("event", &calEventIn);
  RawAnitaHeader* headerIn = NULL;
  headChain->SetBranchAddress("header", &headerIn);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");



  
  //*************************************************************************
  // Set up output
  //*************************************************************************  
  
  TString outFileName = TString::Format("calEventFileWaisPulsesSwappedPolarizations%d.root", firstRun);
  TFile* calEventOutFile = new TFile(outFileName, "recreate");
  TTree* calEventOutTree = new TTree("eventTree", "Tree of Anita Events");
  CalibratedAnitaEvent* calEventOut = NULL;
  calEventOutTree->Branch("event", &calEventOut);

  outFileName = TString::Format("headFileWaisPulsesSwappedPolarizations%d.root", firstRun);
  TFile* headOutFile = new TFile(outFileName, "recreate");
  TTree* headOutTree = new TTree("headTree", "Tree of Anita Headers");
  RawAnitaHeader* headerOut = NULL;

  Double_t waisEventHeading;
  headOutTree->Branch("heading", &waisEventHeading);
  headOutTree->Branch("header", &headerOut);  

  //*************************************************************************
  // Loop over events
  //*************************************************************************  
  
  const Long64_t maxEntries = headChain->GetEntries();
  ProgressBar p(maxEntries);
  for(Long64_t entry=0; entry<maxEntries; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(headerIn->realTime);

    //*************************************************************************
    // WAIS pulse selection based on trigger time
    //*************************************************************************  
    if((headerIn->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      UInt_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      UInt_t triggerTimeNs = headerIn->triggerTimeNs;

      const Double_t maxDeltaTriggerTimeNs = 1200;      
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){    

	// Finally get the big slow thing
	calEventChain->GetEntry(entry);      

	//*************************************************************************
	// Copy header and calibrated event
	//*************************************************************************  
	headerOut = (RawAnitaHeader*) headerIn->Clone();	
	calEventOut = (CalibratedAnitaEvent*) calEventIn->Clone();

	//*************************************************************************
	// Swap event data between V and H channels
	//*************************************************************************  
	for(Int_t polInd=0; polInd<AnitaPol::kNotAPol; polInd++){
	  AnitaPol::AnitaPol_t inputPol = (AnitaPol::AnitaPol_t) polInd;
	  AnitaPol::AnitaPol_t outputPol;
	  if(inputPol==AnitaPol::kHorizontal) {
	    outputPol=AnitaPol::kVertical;
	  }
	  else if(inputPol==AnitaPol::kVertical) {
	    outputPol=AnitaPol::kHorizontal;
	  }
	  else{
	    std::cerr << "??????" << std::endl;
	  }
	  
	  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	    Int_t inputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, inputPol);
	    Int_t outputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, outputPol);

	    Int_t surfIn, chanIn, surfOut, chanOut;
	    AnitaGeomTool::getSurfChanFromChanIndex(inputIndex, surfIn, chanIn);
	    AnitaGeomTool::getSurfChanFromChanIndex(outputIndex, surfOut, chanOut);

	    if(surfIn != surfOut){
	      std::cerr << "Now what am I supposed to do?????" << std::endl;
	    }
	    

	    for(int samp=0; samp<NUM_SAMP; samp++){
	      calEventOut->data[outputIndex][samp] = calEventIn->data[inputIndex][samp];
	    }
	    calEventOut->xMax[outputIndex] = calEventIn->xMax[inputIndex];
	    calEventOut->xMin[outputIndex] = calEventIn->xMin[inputIndex];
	    calEventOut->mean[outputIndex] = calEventIn->mean[inputIndex];
	    calEventOut->rms[outputIndex] = calEventIn->rms[inputIndex];
	  }
	}

	//*************************************************************************
	// Swap header data between V and H branch elements
	//*************************************************************************

	headerOut->l1TrigMask = headerIn->l1TrigMaskH;
	headerOut->l1TrigMaskH = headerIn->l1TrigMask;

	headerOut->phiTrigMask = headerIn->phiTrigMaskH;
	headerOut->phiTrigMaskH = headerIn->phiTrigMask;

	// The lowest bit of the prioritizer is the polarization flag.
	// Because it's the lowest bit I can just add or subtract 1 as required.
	headerOut->prioritizerStuff = headerIn->prioritizerStuff;
	Int_t lowestBit = (headerOut->prioritizerStuff & 1);
	headerOut->prioritizerStuff -= lowestBit; // remove lowest bit
	Int_t newLowestBit = 1 - lowestBit;
	headerOut->prioritizerStuff += newLowestBit; // add new lowest bit

	headerOut->l3TrigPattern = headerIn->l3TrigPatternH;
	headerOut->l3TrigPatternH = headerIn->l3TrigPattern;	


	//*************************************************************************
	// Fill new trees
	//*************************************************************************
	

	waisEventHeading = pat->heading;
	
	headOutTree->Fill();
	calEventOutTree->Fill();
      }
    }
  
    p++;
  }
  calEventOutFile->Write();
  calEventOutFile->Close();  


  headOutFile->Write();
  headOutFile->Close();  

  return 0;
}
