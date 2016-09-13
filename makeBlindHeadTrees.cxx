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


std::vector<std::pair<UInt_t, Int_t> > overwrittenEventInfo; /// pair.first is eventNumber of event to overwrite, pair.second is entry in fakeEventTree to overwrite it with.
TFile* fFakeHeadFile = NULL;
TFile* fFakeEventFile = NULL;
TTree* fFakeHeadTree = NULL;
TTree* fFakeEventTree = NULL;
RawAnitaHeader* fFakeHeader = NULL;
void loadBlindTrees();
Int_t isEventToOverwrite(UInt_t eventNumber);

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


  TString outFileName = TString::Format("blindHeadFileV2_%d.root", firstRun);
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

    Int_t fakeTreeEntry = isEventToOverwrite(headerIn->eventNumber);
    if(fakeTreeEntry >= 0){
      // std::cout << "I'm a fake event!" << headerIn->eventNumber << "\t" << fakeTreeEntry << std::endl;
      // std::cout << fFakeHeadTree << "\t" << fFakeHeader << "\t" << fFakeHeadTree->GetEntry(fakeTreeEntry) << std::endl;
      // std::cout << "the fake info " << fFakeHeader->eventNumber << "\t" << std::endl;
      headerOut = (RawAnitaHeader*) fFakeHeader->Clone();

      // get the ones that obviously stand out in magic display
      headerOut->eventNumber = headerIn->eventNumber;
      headerOut->run = headerIn->run;
      headerOut->trigNum = headerIn->trigNum;
      // std::cout << headerOut->trigNum << std::endl;

      // int getTurfEventNumber()
      // { return (turfEventId&0xfffff);} ///< Returns the event number portion of the TURF event id.
      // headerOut->turfEventId &= ~(0xfffff); // set bits to zero
      // headerOut->turfEventId |= headerIn->getTurfEventNumber(); // set bits to input header
      headerOut->turfEventId = headerIn->turfEventId;
    }
    else{
      headerOut = (RawAnitaHeader*) headerIn->Clone();
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

  // zero internal pointers so can check we find everything.
  fFakeHeadFile = NULL;
  fFakeEventFile = NULL;
  fFakeHeadTree = NULL;
  fFakeEventTree = NULL;
  fFakeHeader = NULL;


  char calibDir[FILENAME_MAX] = ".";
  char fileName[FILENAME_MAX];
  // char *calibEnv=getenv("ANITA_CALIB_DIR");
  // if(!calibEnv) {
  //   char *utilEnv=getenv("ANITA_UTIL_INSTALL_DIR");
  //   if(!utilEnv){
  //     sprintf(calibDir,"calib");
  //   }
  //   else{
  //     sprintf(calibDir,"%s/share/anitaCalib",utilEnv);
  //   }
  // }
  // else {
  //   strncpy(calibDir,calibEnv,FILENAME_MAX);
  // }

  // these are the fake events, that will be inserted in place of some min bias events
  sprintf(fileName, "%s/fakeEventFile.root", calibDir);
  fFakeEventFile = TFile::Open(fileName);
  if(fFakeEventFile){
    fFakeEventTree = (TTree*) fFakeEventFile->Get("eventTree");
  }
  // the header data won't actually get used
  sprintf(fileName, "%s/fakeHeadFile.root", calibDir);
  fFakeHeadFile = TFile::Open(fileName);
  if(fFakeHeadFile){
    fFakeHeadTree = (TTree*) fFakeHeadFile->Get("headTree");
  }


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


  // whinge if you can't find the data
  if(fFakeEventFile && fFakeHeadTree && fFakeHeadFile && fFakeEventTree){
    // std::cerr << "fgInstance = " << fgInstance << ", but this = " << this << std::endl;
    // fFakeEventTree->SetBranchAddress("event", &fFakeEvent);
    fFakeHeadTree->SetBranchAddress("header", &fFakeHeader);
  }
  else{
    std::cerr << "Warning in " << __FILE__ << std::endl;
    std::cerr << "Unable to find files for blinding" << std::endl;
    std::cerr << "fFakeHeadFile = " << fFakeHeadFile << std::endl;
    std::cerr << "fFakeHeadTree = " << fFakeHeadTree << std::endl;
    std::cerr << "fFakeEventFile = " << fFakeEventFile << std::endl;
    std::cerr << "fFakeEventTree = " << fFakeEventTree << std::endl;
  }


}
