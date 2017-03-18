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

int main(int argc, char* argv[]){

  // Runs near WAIS divide
  // const Int_t firstRun = 331;
  // const Int_t lastRun = 354;
  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
  }
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);

  //*************************************************************************
  // Set up input
  //*************************************************************************

  TChain* calEventChain = new TChain("eventTree");
  TChain* headChain = new TChain("headTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    // calEventChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%dOfflineMask.root", run, run);
    headChain->Add(fileName);

  }
  // CalibratedAnitaEvent* calEventIn = NULL;
  // calEventChain->SetBranchAddress("event", &calEventIn);
  RawAnitaHeader* headerIn = NULL;
  headChain->SetBranchAddress("header", &headerIn);


  const int nForTree = 50;
  const UInt_t myPulseEventNumbers[nForTree] = {55602207, 55869718, 55958284, 56017375, 56130483,
						56210753, 56284269, 56355124, 56445987, 56501910,
						56583263, 56697820, 56796435, 56949871, 57094644,
						57209704, 57322092, 57426865, 57519399, 57619903,
						57733132, 57871441, 57980612, 58069692, 58165235,
						58307923, 58439425, 58564882, 58663996, 58728585,
						58794183, 58856764, 58940426, 59111545, 59225479,
						59310321, 59370203, 59558509, 59704990, 59839453,
						60098705, 60258494, 60391416, 60486281, 60563295,
						60630871, 60699867, 60782643, 60917975, 61252049};

  // {55644627, 55870268, 55958489, 56017503, 56130651, 56211033, 56284432, 56355351, 56446163, 56502116, 56583779, 56698192, 56797472, 56950509, 57096157, 57210019, 57322406, 57427501, 57519777, 57620137, 57733483, 57872647, 57981968, 58069949, 58165536, 58308985, 58440526, 58565001, 58664406, 58728758, 58794519, 58857081, 58941263, 59113057, 59225917, 59310637, 59370400, 59558873, 59705951, 59839844, 60100003, 60259393, 60391659, 60486535, 60563477, 60631066, 60700215, 60782984, 60918215, 61273735};



  //*************************************************************************
  // Set up output
  //*************************************************************************

  // TString outFileName = "fakeEventFile.root";
  // TFile* usefulEventOutFile = new TFile(outFileName, "recreate");
  // TTree* usefulEventOutTree = new TTree("eventTree", "Tree of Anita Events");
  // UsefulAnitaEvent* usefulEventOut = NULL;
  // usefulEventOutTree->Branch("event", &usefulEventOut);

  TString outFileName = "fakeHeadFile.root";
  TFile* headOutFile = new TFile(outFileName, "recreate");
  TTree* headOutTree = new TTree("headTree", "Tree of Anita Headers");
  RawAnitaHeader* headerOut = NULL;

  Double_t waisEventHeading;
  headOutTree->Branch("heading", &waisEventHeading);
  headOutTree->Branch("header", &headerOut);


  //*************************************************************************
  // Loop over pulse indices
  //*************************************************************************

  headChain->BuildIndex("eventNumber");

  const Long64_t maxEntries = nForTree;
  ProgressBar p(maxEntries);
  for(Long64_t pulseInd=0; pulseInd<nForTree; pulseInd++){
    Long64_t entry = headChain->GetEntryNumberWithIndex(myPulseEventNumbers[pulseInd]);

    // std::cout << entry << std::endl;
    if(entry >= 0){
      headChain->GetEntry(entry);
      calEventChain->GetEntry(entry);

      //*************************************************************************
      // Copy header and calibrated event
      //*************************************************************************
      // headerOut = (RawAnitaHeader*) headerIn->Clone();
      headerOut = headerIn; // copy pointer



      // usefulEventOut = (CalibratedAnitaEvent*) calEventIn->Clone();
      // UsefulAnitaEvent* usefulEventTemp = new UsefulAnitaEvent(calEventIn);
      // usefulEventOut = new UsefulAnitaEvent(calEventIn);

      //*************************************************************************
      // Swap event data between V and H channels
      //*************************************************************************

      // for(Int_t polInd=0; polInd<AnitaPol::kNotAPol; polInd++){
      // 	AnitaPol::AnitaPol_t inputPol = (AnitaPol::AnitaPol_t) polInd;
      // 	AnitaPol::AnitaPol_t outputPol;
      // 	if(inputPol==AnitaPol::kHorizontal) {
      // 	  outputPol=AnitaPol::kVertical;
      // 	}
      // 	else if(inputPol==AnitaPol::kVertical) {
      // 	  outputPol=AnitaPol::kHorizontal;
      // 	}
      // 	else{
      // 	  std::cerr << "??????" << std::endl;
      // 	}

      // 	for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      // 	  Int_t inputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, inputPol);
      // 	  Int_t outputIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, outputPol);

      // 	  Int_t surfIn, chanIn, surfOut, chanOut;
      // 	  AnitaGeomTool::getSurfChanFromChanIndex(inputIndex, surfIn, chanIn);
      // 	  AnitaGeomTool::getSurfChanFromChanIndex(outputIndex, surfOut, chanOut);

      // 	  if(surfIn != surfOut){
      // 	    std::cerr << "Now what am I supposed to do?????" << std::endl;
      // 	  }

      // 	  for(int samp=0; samp<NUM_SAMP; samp++){
      // 	    usefulEventOut->data[outputIndex][samp] = usefulEventTemp->data[inputIndex][samp];
      // 	  }
      // 	  usefulEventOut->xMax[outputIndex] = usefulEventTemp->xMax[inputIndex];
      // 	  usefulEventOut->xMin[outputIndex] = usefulEventTemp->xMin[inputIndex];
      // 	  usefulEventOut->mean[outputIndex] = usefulEventTemp->mean[inputIndex];
      // 	  usefulEventOut->rms[outputIndex] = usefulEventTemp->rms[inputIndex];

      // 	  usefulEventOut->fNumPoints[outputIndex] = usefulEventTemp->fNumPoints[inputIndex];

      // 	  for(int samp=0; samp < NUM_SAMP; samp++){
      // 	    usefulEventOut->fVolts[outputIndex][samp] = usefulEventTemp->fVolts[inputIndex][samp];
      // 	    usefulEventOut->fTimes[outputIndex][samp] = usefulEventTemp->fTimes[inputIndex][samp];
      // 	  }

      // 	  // if input is ALFA
      // 	  if(inputIndex == (11*NUM_CHAN + 5)){
      // 	    // std::cout << usefulEventOut->fVolts[outputIndex][0] << "\t"
      // 	    // 	      << usefulEventOut->fVolts[outputIndex][1] << "\t"
      // 	    // 	      << usefulEventOut->fVolts[outputIndex][2] << "\t"
      // 	    // 	      << std::endl;

      // 	    std::complex<double>* theFFT = FancyFFTs::doFFT(usefulEventOut->fNumPoints[outputIndex],
      // 							    &usefulEventOut->fVolts[outputIndex][0],
      // 							    true);

      // 	    const int nf = FancyFFTs::getNumFreqs(usefulEventOut->fNumPoints[outputIndex]);
      // 	    double deltaF = 1e3/((1./2.6)*usefulEventOut->fNumPoints[outputIndex]);

      // 	    for(int i=0; i < nf; i++){
      // 	      double freq = deltaF*i;
      // 	      if(freq >= 700){
      // 		theFFT[i].real(0);
      // 		theFFT[i].imag(0);
      // 	      }
      // 	    }

      // 	    FancyFFTs::doInvFFT(usefulEventOut->fNumPoints[outputIndex],
      // 				theFFT,
      // 				&usefulEventOut->fVolts[outputIndex][0],
      // 				true);
      // 	    delete [] theFFT;

      // 	    // std::cout << usefulEventOut->fVolts[outputIndex][0] << "\t"
      // 	    // 	      << usefulEventOut->fVolts[outputIndex][1] << "\t"
      // 	    // 	      << usefulEventOut->fVolts[outputIndex][2] << "\t"
      // 	    // 	      << std::endl;
      // 	  }
      // 	}
      // }

      //*************************************************************************
      // Swap header data between V and H branch elements
      //*************************************************************************

      // headerOut->l1TrigMask = headerIn->l1TrigMaskH;
      // headerOut->l1TrigMaskH = headerIn->l1TrigMask;

      // headerOut->phiTrigMask = headerIn->phiTrigMaskH;
      // headerOut->phiTrigMaskH = headerIn->phiTrigMask;

      // The lowest bit of the prioritizer is the polarization flag.
      // Because it's the lowest bit I can just add or subtract 1 as required.
      // headerOut->prioritizerStuff = headerIn->prioritizerStuff;
      // Int_t lowestBit = (headerOut->prioritizerStuff & 1);
      // headerOut->prioritizerStuff -= lowestBit; // remove lowest bit
      // Int_t newLowestBit = 1 - lowestBit;
      // headerOut->prioritizerStuff += newLowestBit; // add new lowest bit

      // headerOut->l3TrigPattern = headerIn->l3TrigPatternH;
      // headerOut->l3TrigPatternH = headerIn->l3TrigPattern;

      //*************************************************************************
      // Fill new trees
      //*************************************************************************

      headOutTree->Fill();
      // usefulEventOutTree->Fill();

      // delete usefulEventOut;
      // usefulEventOut = NULL;

      // delete usefulEventTemp;
      // waisEventHeading = pat->heading;
    }

    p.inc(pulseInd, nForTree);
  }
  // usefulEventOutFile->Write();
  // usefulEventOutFile->Close();


  headOutFile->Write();
  headOutFile->Close();

  return 0;
}
