// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Selects some minimum bias events (not in the decimated data set) from a set of time bins spanning the flight, and writes them to a file. Pretty simple really.
*************************************************************************************************************** */

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TApplication.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"
#include "AnitaVersion.h"

#include "AnitaEventSummary.h"
#include "OutputConvention.h"
#include "ProgressBar.h"
#include "TGraphAntarctica.h"
#include "FancyFFTs.h"
#include "RootTools.h"


int main(int argc, char* argv[]){

  AnitaVersion::set(3);

  TApplication* theApp = new TApplication(argv[0], &argc, argv);

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
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* decimated = new TChain("headTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run < 257 || run > 263){
      TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/timedHeadFile%dOfflineMask.root", run, run);
      headChain->Add(fileName);

      fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
      gpsChain->Add(fileName);

      fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
      decimated->Add(fileName);
    }
  }


  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  headChain->BuildIndex("realTime");
  decimated->BuildIndex("eventNumber");

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  const int numPulses=50;

  std::ofstream outFile("anita3OverwrittenEventInfo.txt", std::ofstream::out);
  outFile << "eventNumber\tfakeTreeEntry" << std::endl;

  TFile* fReco = OutputConvention::getFile("reconstruction_*");
  TTree* tReco = (TTree*) fReco->Get("eventSummaryTree");
  AnitaEventSummary* summary = NULL;
  tReco->SetBranchAddress("eventSummary", &summary);


  // Now for version 3 of the blinding, we will divide the flight into segments of time
  // rather than segments of eventNumber...
  headChain->GetEntry(0);
  UInt_t firstRealTime = header->realTime;
  headChain->GetEntry(headChain->GetEntries()-1);
  UInt_t lastRealTime = header->realTime;
  // std::cout << firstRealTime << "\t" << lastRealTime << "\t" << lastRealTime - firstRealTime << std::endl;

  UInt_t seed = 29348756; // mashed keyboard with hands
  // UInt_t seed = 13986513; // mashed keyboard with hands
  // UInt_t seed = 0;
  TRandom3 rnd(seed);

  // don't print this number on the final go...
  const int N = rnd.Uniform(10, 15);
  // const int N = 10;

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


  // Long64_t nEntries = headChain->GetEntries();

  // now loop over number of selected blinded pulses, and place them

  TGraphAntarctica* grBlindRecoPosition = new TGraphAntarctica();
  TGraphAntarctica* grAnitaPat = new TGraphAntarctica();
  std::vector<TGraphAntarctica*> grConnectors;
  Double_t eventBearing = -1;
  // ProgressBar p(N);
  for(Long64_t i=0; i < N; i++){



    // get the fake event
    Int_t fakeTreeEntry = fakeTreeEntries.at(selectedEvents);
    tReco->GetEntry(fakeTreeEntry);

    Int_t isDec = 1;
    Int_t isMinBias = 0;

    int numWhile = 0;

    bool onContinent = false;
    bool southFacingEvent = false;

    Double_t sourceLon, sourceLat, sourceAltitude;

    const int crazyNumber = 1000000;
    while(isDec > -1 || isMinBias <= 0 || !onContinent || !southFacingEvent){

      // Pick any event from within the time window
      UInt_t randomTime = rnd.Uniform(firstRealTime, lastRealTime);
      Int_t entry = headChain->GetEntryNumberWithIndex(randomTime);

      headChain->GetEntry(entry);
      isDec = decimated->GetEntryNumberWithIndex(header->eventNumber);
      isMinBias = header->getTriggerBitSoftExt();



      // Now check if reconstructs to the continent
      gpsChain->GetEntry(entry);
      // std::cout << gpsBytes << std::endl;
      UsefulAdu5Pat usefulPat(pat);


      Double_t phiWave = summary->peak[AnitaPol::kVertical][0].phi*TMath::DegToRad();
      Double_t thetaWave = summary->peak[AnitaPol::kVertical][0].theta*TMath::DegToRad();
      int retVal = usefulPat.getSourceLonAndLatAtAlt(phiWave, -thetaWave, sourceLon, sourceLat, sourceAltitude);

      if(retVal!=1){
	onContinent = false;
      }
      else{
	// skip events that don't reconstruct to continent
	onContinent = RampdemReader::isOnContinent(sourceLon, sourceLat);
      }

      eventBearing = RootTools::getDeltaAngleDeg(pat->heading, summary->peak[AnitaPol::kVertical][0].phi);
      southFacingEvent = TMath::Abs(eventBearing) > 135;
      std::cout << eventBearing << "\t" << pat->heading << "\t" << summary->peak[AnitaPol::kVertical][0].phi << std::endl;


      // std::cout << "In while loop " << isDec << "\t" << isMinBias << "\t" << onContinent << std::endl;

      if(numWhile >= crazyNumber){
	std::cerr << "Something seems wrong here. I'm quitting!" << std::endl;
	return 1;
      }
      numWhile++;

    }

    grBlindRecoPosition->SetPoint(grBlindRecoPosition->GetN(), sourceLon, sourceLat);
    grAnitaPat->SetPoint(grAnitaPat->GetN(), pat->longitude, pat->latitude);

    TGraphAntarctica* grConnector = new TGraphAntarctica();
    grConnector->SetPoint(grConnector->GetN(), pat->longitude, pat->latitude);
    grConnector->SetPoint(grConnector->GetN(), sourceLon, sourceLat);
    grConnectors.push_back(grConnector);

    UsefulAdu5Pat usefulPat(pat);

    Double_t sourceAlt = RampdemReader::SurfaceAboveGeoid(sourceLon, sourceLat);
    Double_t distKm = 1e-3*usefulPat.getDistanceFromSource(sourceLat, sourceLon, sourceAlt);

    std::cout << "Inserted event " << i << ":" << std::endl;
    std::cout << "ANITA at " << pat->longitude << "\t" << pat->latitude << "\t" << 1e-3*pat->altitude << std::endl;
    std::cout << "Reconstructed position at " << sourceLon << "\t" << sourceLat << "\t" << 1e-3*sourceAlt << std::endl;
    std::cout << "They are separated by " << distKm << " km"  << std::endl;
    std::cout << "Event bearing = " << eventBearing << "\t" << pat->heading << "\t" << summary->peak[AnitaPol::kVertical][0].phi << std::endl;

    // write event number to file
    outFile << header->eventNumber << "\t" << fakeTreeEntry << std::endl;

    selectedEvents++;
    // p.inc(i, N);
    std::cout << i << "\t" << N << std::endl;
  }

  outFile.close();

  grBlindRecoPosition->SetMarkerStyle(8);
  grBlindRecoPosition->SetMarkerColor(kRed);
  grBlindRecoPosition->Draw("ap");

  const int numPointsIWant = 5000;
  const int selectEvery = (gpsChain->GetEntries())/numPointsIWant;
  TString cutString = "(Entry$ % " + TString::Format("%d)==0", selectEvery);

  TGraphAntarctica* grFlightPath = new TGraphAntarctica(gpsChain, "longitude", "latitude", TCut(cutString.Data()));
  grFlightPath->SetLineColor(kGreen);
  grFlightPath->Draw("lsame");

  grAnitaPat->SetMarkerStyle(8);
  grAnitaPat->SetMarkerColor(kGreen);
  grAnitaPat->Draw("psame");

  TLegend* l = new TLegend(0.75, 0, 1, 0.25);
  l->AddEntry(grFlightPath, "ANITA-3 Flight Path", "l");
  l->AddEntry(grAnitaPat, "ANITA's position for inserted events", "p");
  l->AddEntry(grBlindRecoPosition, "Reconstructed position of inserted events", "p");
  l->AddEntry(grConnectors.at(0), "Reconstructed position to ANITA", "l");


  for(UInt_t i=0; i < grConnectors.size(); i++){
    grConnectors.at(i)->Draw("lsame");
    grConnectors.at(i)->SetLineColor(kMagenta);
  }
  l->Draw();

  theApp->Run();

  return 0;

}
