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
#include "FancyFFTs.h"


TFile* fakeEventFile = NULL;
TTree* fakeEventTree = NULL;
UsefulAnitaEvent* fakeEvent = NULL;

// pair.first is eventNumber of event to overwrite
// pair.second is entry in fakeEventTree to overwrite it with.
std::vector<std::pair<UInt_t, Int_t> > overwrittenEventInfo;

void loadBlindTrees();
Int_t isEventToOverwrite(UInt_t eventNumber);

Int_t blindingVersion = 3; // since finishing thesis

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

  loadBlindTrees();

  TChain* headChain = new TChain("headTree");
  TChain* fakeChain = new TChain("headTree");
  for(Int_t run=331; run <= 354; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%dOfflineMask.root", run, run);
    fakeChain->Add(fileName);
  }
  fakeChain->BuildIndex("eventNumber");
  RawAnitaHeader* fakeHeader = NULL;
  fakeChain->SetBranchAddress("header", &fakeHeader);


  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%dOfflineMask.root", run, run);
    headChain->Add(fileName);

  }
  RawAnitaHeader* headerIn = NULL;
  headChain->SetBranchAddress("header", &headerIn);

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header file for run " << firstRun << ". Giving up." << std::endl;
    return 1;
  }

  //*************************************************************************
  // Set up output
  //*************************************************************************


  TString outFileName = TString::Format("blindHeadFileV%d_%d.root", blindingVersion, firstRun);
  TFile* headOutFile = new TFile(outFileName, "recreate");
  TTree* headOutTree = new TTree("headTree", "Tree of Anita Headers");
  RawAnitaHeader* headerOut = NULL;
  headOutTree->Branch("header", &headerOut);
  //*************************************************************************
  // Loop over pulse indices
  //*************************************************************************

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry=startEntry; entry<maxEntry; entry++){

    headChain->GetEntry(entry);

    //*************************************************************************
    // Copy header and calibrated event
    //*************************************************************************


    headerOut = headerIn;

    Int_t fakeTreeEntry = isEventToOverwrite(headerIn->eventNumber);
    if(fakeTreeEntry >= 0){

      fakeEventTree->GetEntry(fakeTreeEntry);
      fakeChain->GetEntryWithIndex(fakeEvent->eventNumber);

      headerOut->l1TrigMask = fakeHeader->l1TrigMaskH;
      headerOut->l1TrigMaskH = fakeHeader->l1TrigMask;
      headerOut->phiTrigMask = fakeHeader->phiTrigMaskH;
      headerOut->phiTrigMaskH = fakeHeader->phiTrigMask;
      headerOut->l1TrigMaskOffline = fakeHeader->l1TrigMaskHOffline;
      headerOut->l1TrigMaskHOffline = fakeHeader->l1TrigMaskOffline;
      headerOut->phiTrigMaskOffline = fakeHeader->phiTrigMaskHOffline;
      headerOut->phiTrigMaskHOffline = fakeHeader->phiTrigMaskOffline;


      headerOut->l3TrigPattern = fakeHeader->l3TrigPatternH;
      headerOut->l3TrigPatternH = fakeHeader->l3TrigPattern;

      // looks like TObject::Clone doesn't properly copy UChar_t (maybe Char_t too?)
      // so do this manaully here
      headerOut->priority = fakeHeader->priority;
      headerOut->turfUpperWord = fakeHeader->turfUpperWord;
      headerOut->otherFlag = fakeHeader->otherFlag;
      headerOut->errorFlag = fakeHeader->errorFlag;
      headerOut->surfSlipFlag = fakeHeader->surfSlipFlag = fakeHeader->surfSlipFlag = fakeHeader->surfSlipFlag;
      headerOut->nadirAntTrigMask = fakeHeader->nadirAntTrigMask;
      headerOut->peakThetaBin = fakeHeader->peakThetaBin;
      for(int i=0; i < 2; i++){
	headerOut->reserved[i] = fakeHeader->reserved[i];
      }
      headerOut->trigType = fakeHeader->trigType;
      headerOut->l3Type1Count = fakeHeader->l3Type1Count;
      headerOut->bufferDepth = fakeHeader->bufferDepth;
      headerOut->turfioReserved = fakeHeader->turfioReserved;
      headerOut->nadirL1TrigPattern = fakeHeader->nadirL1TrigPattern;
      headerOut->nadirL2TrigPattern = fakeHeader->nadirL2TrigPattern;



      std::cout << headerOut->eventNumber << "\t" << headerOut->trigNum << std::endl;
      std::cout << (fakeHeader->errorFlag & 0x1) << "\t" << (headerOut->errorFlag & 0x1) << std::endl;
      std::cout << (fakeHeader->errorFlag & 0x2) << "\t" << (headerOut->errorFlag & 0x2) << std::endl;
      std::cout << (fakeHeader->errorFlag & 0x4) << "\t" << (headerOut->errorFlag & 0x4) << std::endl;
      std::cout << (fakeHeader->errorFlag & 0x8) << "\t" << (headerOut->errorFlag & 0x8) << std::endl;
      std::cout << (fakeHeader->errorFlag & 0x10) << "\t" << (headerOut->errorFlag & 0xf) << std::endl;

      std::cout << (fakeHeader->priority) << "\t" << (headerOut->priority) << std::endl;
      std::cout << (fakeHeader->turfUpperWord) << "\t" << (headerOut->turfUpperWord) << std::endl;
      std::cout << (fakeHeader->otherFlag) << "\t" << (headerOut->otherFlag) << std::endl;
      std::cout << (fakeHeader->surfSlipFlag) << "\t" << (headerOut->surfSlipFlag) << std::endl;

      fakeEvent = NULL;

    }

    headOutTree->Fill();

    p.inc(entry, maxEntry);
  }
  headOutFile->Write();
  headOutFile->Close();

  return 0;
}


Int_t isEventToOverwrite(UInt_t eventNumber){

  Int_t fakeTreeEntry = -1;
  for(UInt_t i=0; i <overwrittenEventInfo.size(); i++){
    if(overwrittenEventInfo.at(i).first==eventNumber){
      fakeTreeEntry = overwrittenEventInfo.at(i).second;
      break;
    }
  }
  return fakeTreeEntry;
}


void loadBlindTrees() {

  char calibDir[FILENAME_MAX] = ".";
  char fileName[FILENAME_MAX];



  // these are the min bias event numbers to be overwritten, with the entry in the fakeEventTree
  // that is used to overwrite the event
  sprintf(fileName,"%s/anita3OverwrittenEventInfo.txt",calibDir);
  std::ifstream overwrittenEventInfoFile(fileName);
  char firstLine[180];
  overwrittenEventInfoFile.getline(firstLine,179);
  UInt_t overwrittenEventNumber;
  Int_t fakeTreeEntry;
  while(overwrittenEventInfoFile >> overwrittenEventNumber >> fakeTreeEntry){
    overwrittenEventInfo.push_back(std::pair<UInt_t, Int_t>(overwrittenEventNumber, fakeTreeEntry));
    // std::cout << overwrittenEventInfo.at(overwrittenEventInfo.size()-1).first << "\t" << overwrittenEventInfo.at(overwrittenEventInfo.size()-1).second << std::endl;
  }
  if(overwrittenEventInfo.size()==0){
    std::cerr << "Warning in " << __FILE__ << std::endl;
    std::cerr << "Unable to find overwrittenEventInfo" << std::endl;
  }


  fakeEventFile = TFile::Open("fakeEventFile.root");
  fakeEventTree = (TTree*) fakeEventFile->Get("eventTree");
  fakeEventTree->SetBranchAddress("event", &fakeEvent);

}
