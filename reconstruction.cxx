// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Reconstruct entire data set.
********************************************************************************************************* */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"
#include "AnitaDataSet.h"

int main(int argc, char *argv[]){


  const Int_t firstRun = 331;
  const Int_t lastRun = 354;

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventChain = new TChain("eventTree");

  const char* anitaInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  if(anitaInstallDir==NULL){
    std::cerr << "you fool" << std::endl;
    exit(1);
  }
  TString headerFileName = TString::Format("%s/share/anitaCalib/fakeHeadFile.root", anitaInstallDir);
  headChain->Add(headerFileName);
  TString eventFileName = TString::Format("%s/share/anitaCalib/fakeHeadFile.root", anitaInstallDir);
  headChain->Add(eventFileName);

  for(Int_t run=firstRun; run<=lastRun; run++){

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
  }

  gpsChain->BuildIndex("eventNumber");

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  UsefulAnitaEvent* usefulEvent = NULL;
  eventChain->SetBranchAddress("event", &usefulEvent);

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite And 200MHz Notch Notch",
  					260-26, 260+26);
  CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch",
  					370-26, 370+26);
  CrossCorrelator::SimpleNotch notch400("n400Notch", "400 MHz Satellite Notch",
  					 400-10, 410);
  CrossCorrelator::SimpleNotch notch762("n762Notch", "762MHz Satellite Notch (one bin wide)",
  					762-8, 762+8);
  CrossCorrelator::SimpleNotch notch200("n200Notch", "200 MHz high pass band",
  					0, 200);
  CrossCorrelator::SimpleNotch notch1200("n1200Notch", "1200 MHz low pass band",
  					 1200, 9999);

  cc->addNotch(notch260);
  cc->addNotch(notch370);
  cc->addNotch(notch400);
  cc->addNotch(notch762);
  cc->addNotch(notch200);
  cc->addNotch(notch1200);

  const Int_t myNumPeaksCoarse = 5;
  const Int_t myNumPeaksFine = 5;
  const Int_t coherentDeltaPhi = 0;

  TTree* eventSummaryTree = new TTree("eventSummaryTree", "eventSummaryTree");
  // AnitaEventSummary* eventSummary = new AnitaEventSummary();
  AnitaEventSummary* eventSummary = NULL; //new AnitaEventSummary();
  eventSummaryTree->Branch("eventSummary", &eventSummary);


  Long64_t nEntries = headChain->GetEntries();
  Long64_t startEntry = 0;
  Long64_t maxEntry = headChain->GetEntries();

  std::cout << "Processing " << maxEntry-startEntry << " of " << nEntries << " entries." << std::endl;
  std::cout << "Starting at entry " << startEntry << " up to  " << maxEntry << " entry." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);
    eventChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->eventNumber);

    UsefulAdu5Pat usefulPat(pat);
    cc->reconstructEvent(usefulEvent, myNumPeaksCoarse, myNumPeaksFine);

    eventSummary = new AnitaEventSummary(header, &usefulPat);
    // std::cout << eventSummary->sun.theta << "\t" << eventSummary->sun.phi << std::endl;

    Double_t minY = 0;

    for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

      for(Int_t peakInd=0; peakInd < myNumPeaksFine; peakInd++){
	cc->getFinePeakInfo(pol, peakInd,
			    eventSummary->peak[pol][peakInd].value,
			    eventSummary->peak[pol][peakInd].phi,
			    eventSummary->peak[pol][peakInd].theta);

	TGraph* grZ0 = cc->makeUpsampledCoherentlySummedWaveform(pol,
								 eventSummary->peak[pol][peakInd].phi,
								 eventSummary->peak[pol][peakInd].theta,
								 coherentDeltaPhi,
								 // eventSummary->peak[pol][peakInd].snr);
								 eventSummary->coherent[pol][peakInd].snr);

	if(grZ0!=NULL){
	  TGraph* grZ0Hilbert = FFTtools::getHilbertEnvelope(grZ0);

	  RootTools::getMaxMin(grZ0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);

	  delete grZ0;
	  delete grZ0Hilbert;
	}
      }
    }

    eventSummary->flags.isGood = 1;
    eventSummary->flags.isPayloadBlast = 0; //!< To be determined.
    eventSummary->flags.nadirFlag = 0; //!< Not sure I will use this.
    eventSummary->flags.strongCWFlag = 0; //!< Not sure I will use this.
    eventSummary->flags.isVarner = 0; //!< Not sure I will use this.
    eventSummary->flags.isVarner2 = 0; //!< Not sure I will use this.
    eventSummary->flags.pulser = AnitaEventSummary::EventFlags::NONE; //!< Not yet.

    delete usefulEvent;

    eventSummaryTree->Fill();
    // delete eventSummary;
    p.inc(entry, nEntries);
  }

  // saves time later
  eventSummaryTree->BuildIndex("eventNumber");

  outFile->Write();
  outFile->Close();

  return 0;
}
