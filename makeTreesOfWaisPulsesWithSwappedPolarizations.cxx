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
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%dOfflineMask.root", run, run);
    headChain->Add(fileName);

  }

  headChain->BuildIndex("eventNumber");

  CalibratedAnitaEvent* calEventIn = NULL;
  calEventChain->SetBranchAddress("event", &calEventIn);
  RawAnitaHeader* headerIn = NULL;
  headChain->SetBranchAddress("header", &headerIn);


  // const int nForTree = 50;
  const int nPerTree = 25;
  const int numPol = 2;
  const UInt_t myPulseEventNumbers[numPol][nPerTree] = {{55602207, 55869718, 55958284, 56017375, 56130483,
							 56210753, 56284269, 56355124, 56445987, 56501910,
							 56583263, 56697820, 56796435, 56949871, 57094644,
							 57209704, 57322092, 57426865, 57519399, 57619903,
							 57733132, 57871441, 57980612, 58069692, 58165235},
							{58307923, 58439425, 58564882, 58663996, 58728585,
							 58794183, 58856764, 58940426, 59111545, 59225479,
							 59310321, 59370203, 59558509, 59704990, 59839453,
							 60098705, 60258494, 60391416, 60486281, 60563295,
							 60630871, 60699867, 60782643, 60917975, 61252049}};


  //*************************************************************************
  // Set up output
  //*************************************************************************

  TString outFileName = "headersAndEventsForBlindingByInsertion.root";
  // TString outFileName = "insertedEvents.root";
  // TString outFileName = "fakeEventFile.root";
  TFile* outFile = new TFile(outFileName, "recreate");

  TString polNames[AnitaPol::kNotAPol] = {"HPol", "VPol"};
  TTree* outEventTrees[AnitaPol::kNotAPol] = {NULL};
  TTree* outHeadTrees[AnitaPol::kNotAPol] = {NULL};

  for(int polIndTree=0; polIndTree < AnitaPol::kNotAPol; polIndTree++){
    outEventTrees[polIndTree] = new TTree(polNames[polIndTree] + "EventTree", "Tree of Anita Events");
    outHeadTrees[polIndTree] = new TTree(polNames[polIndTree] + "HeadTree", "Tree of Anita Headers");

    TTree* usefulEventOutTree = outEventTrees[polIndTree];
    TTree* headOutTree = outHeadTrees[polIndTree];

    UsefulAnitaEvent* usefulEventOut = NULL;
    usefulEventOutTree->Branch("event", &usefulEventOut);


    RawAnitaHeader* headerOut = NULL;
    Double_t waisEventHeading;
    headOutTree->Branch("heading", &waisEventHeading);
    headOutTree->Branch("header", &headerOut);





    //*************************************************************************
    // Loop over pulse indices
    //*************************************************************************

    const Long64_t maxEntries = nPerTree;
    ProgressBar p(maxEntries);
    for(Long64_t pulseInd=0; pulseInd<nPerTree; pulseInd++){
      Long64_t entry = headChain->GetEntryNumberWithIndex(myPulseEventNumbers[polIndTree][pulseInd]);

      // std::cout << entry << std::endl;
      if(entry >= 0){
	headChain->GetEntry(entry);
	calEventChain->GetEntry(entry);

	//*************************************************************************
	// Copy header and calibrated event
	//*************************************************************************
	headerOut = headerIn; // copy pointer

	UsefulAnitaEvent* usefulEventTemp = new UsefulAnitaEvent(calEventIn);
	usefulEventOut = new UsefulAnitaEvent(calEventIn);

	if(polIndTree==AnitaPol::kVertical){

	  // *************************************************************************
	  // Swap event data between V and H channels
	  // *************************************************************************

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
		usefulEventOut->data[outputIndex][samp] = usefulEventTemp->data[inputIndex][samp];
	      }
	      usefulEventOut->xMax[outputIndex] = usefulEventTemp->xMax[inputIndex];
	      usefulEventOut->xMin[outputIndex] = usefulEventTemp->xMin[inputIndex];
	      usefulEventOut->mean[outputIndex] = usefulEventTemp->mean[inputIndex];
	      usefulEventOut->rms[outputIndex] = usefulEventTemp->rms[inputIndex];

	      usefulEventOut->fNumPoints[outputIndex] = usefulEventTemp->fNumPoints[inputIndex];

	      for(int samp=0; samp < NUM_SAMP; samp++){
		usefulEventOut->fVolts[outputIndex][samp] = usefulEventTemp->fVolts[inputIndex][samp];
		usefulEventOut->fTimes[outputIndex][samp] = usefulEventTemp->fTimes[inputIndex][samp];
	      }

	      // if input is ALFA
	      if(inputIndex == (11*NUM_CHAN + 5)){
		// std::cout << usefulEventOut->fVolts[outputIndex][0] << "\t"
		// 	      << usefulEventOut->fVolts[outputIndex][1] << "\t"
		// 	      << usefulEventOut->fVolts[outputIndex][2] << "\t"
		// 	      << std::endl;

		std::complex<double>* theFFT = FancyFFTs::doFFT(usefulEventOut->fNumPoints[outputIndex],
								&usefulEventOut->fVolts[outputIndex][0],
								true);

		const int nf = FancyFFTs::getNumFreqs(usefulEventOut->fNumPoints[outputIndex]);
		double deltaF = 1e3/((1./2.6)*usefulEventOut->fNumPoints[outputIndex]);

		for(int i=0; i < nf; i++){
		  double freq = deltaF*i;
		  if(freq >= 700){
		    theFFT[i].real(0);
		    theFFT[i].imag(0);
		  }
		}

		FancyFFTs::doInvFFT(usefulEventOut->fNumPoints[outputIndex],
				    theFFT,
				    &usefulEventOut->fVolts[outputIndex][0],
				    true);
		delete [] theFFT;

		// std::cout << usefulEventOut->fVolts[outputIndex][0] << "\t"
		// 	      << usefulEventOut->fVolts[outputIndex][1] << "\t"
		// 	      << usefulEventOut->fVolts[outputIndex][2] << "\t"
		// 	      << std::endl;
	      }
	    }
	  }

	  RawAnitaHeader fakeHeader2 = (*headerOut);
	  RawAnitaHeader* fakeHeader = &fakeHeader2;

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

	}


	//*************************************************************************
	// Fill new trees
	//*************************************************************************

	headOutTree->Fill();
	usefulEventOutTree->Fill();

	delete usefulEventTemp;
	delete usefulEventOut;
	usefulEventOut = NULL;

      }

      p.inc(pulseInd, nPerTree);
    }
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
