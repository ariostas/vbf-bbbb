#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

void cleanUpMergedFiles(TString infilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_temp.root",
			  TString outfilename="/afs/cern.ch/work/a/ariostas/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_clean.root") {

  // set up input/output variables and file
  UInt_t nEvents;
  Double_t weight;
  UInt_t nLeptons=0, nLeptons01=0, nLeptons04=0, nJets=0, nBJets=0, nJetsPU=0, nBJetsPU=0, nDiJets=0, matchedBJets=0, matchedLJets=0, matchedHiggs=0;
  
  Double_t bjet1_Pt, bjet1_Eta, bjet1_Phi, bjet1_Mass;
	Double_t bjet2_Pt, bjet2_Eta, bjet2_Phi, bjet2_Mass;
	Double_t bjet3_Pt, bjet3_Eta, bjet3_Phi, bjet3_Mass;
	Double_t bjet4_Pt, bjet4_Eta, bjet4_Phi, bjet4_Mass;
	Double_t jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass;
	Double_t jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass;

  TFile* infile = new TFile(infilename); assert(infile);
  TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("weight",		&weight);
  intree->SetBranchAddress("bjet1_Pt",		&bjet1_Pt);
  intree->SetBranchAddress("bjet1_Eta",		&bjet1_Eta);
  intree->SetBranchAddress("bjet1_Phi",		&bjet1_Phi);
  intree->SetBranchAddress("bjet1_Mass",	&bjet1_Mass);
  intree->SetBranchAddress("bjet2_Pt",		&bjet2_Pt);
  intree->SetBranchAddress("bjet2_Eta",		&bjet2_Eta);
  intree->SetBranchAddress("bjet2_Phi",		&bjet2_Phi);
  intree->SetBranchAddress("bjet2_Mass",	&bjet2_Mass);
  intree->SetBranchAddress("bjet3_Pt",		&bjet3_Pt);
  intree->SetBranchAddress("bjet3_Eta",		&bjet3_Eta);
  intree->SetBranchAddress("bjet3_Phi",		&bjet3_Phi);
  intree->SetBranchAddress("bjet3_Mass",	&bjet3_Mass);
  intree->SetBranchAddress("bjet4_Pt",		&bjet4_Pt);
  intree->SetBranchAddress("bjet4_Eta",		&bjet4_Eta);
  intree->SetBranchAddress("bjet4_Phi",		&bjet4_Phi);
  intree->SetBranchAddress("bjet4_Mass",	&bjet4_Mass);
  intree->SetBranchAddress("jet1_Pt",		&jet1_Pt);
  intree->SetBranchAddress("jet1_Eta",		&jet1_Eta);
  intree->SetBranchAddress("jet1_Phi",		&jet1_Phi);
  intree->SetBranchAddress("jet1_Mass",		&jet1_Mass);
  intree->SetBranchAddress("jet2_Pt",		&jet2_Pt);
  intree->SetBranchAddress("jet2_Eta",		&jet2_Eta);
  intree->SetBranchAddress("jet2_Phi",		&jet2_Phi);
  intree->SetBranchAddress("jet2_Mass",		&jet2_Mass);
  intree->SetBranchAddress("nLeptons",     	&nLeptons);
  intree->SetBranchAddress("nLeptons01",    &nLeptons01);
  intree->SetBranchAddress("nLeptons04",    &nLeptons04);
  intree->SetBranchAddress("nJets",     	&nJets);
  intree->SetBranchAddress("nBJets",     	&nBJets);
  intree->SetBranchAddress("nJetsPU",     	&nJetsPU);
  intree->SetBranchAddress("nBJetsPU",     	&nBJetsPU);
  intree->SetBranchAddress("nDiJets",     	&nDiJets);
  intree->SetBranchAddress("matchedBJets",  &matchedBJets);
  intree->SetBranchAddress("matchedLJets",  &matchedLJets);
  intree->SetBranchAddress("matchedHiggs",  &matchedHiggs);

  TTree* infotree = (TTree*) infile->Get("Info"); assert(infotree);
  infotree->SetBranchAddress("nEvents",      &nEvents);

  Long64_t totalEvents=0;

  for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
    infotree->GetEntry(iEntry);
    totalEvents+=nEvents;
  }

  TFile *outFile = new TFile(outfilename, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",				&weight,    			"weight/D");  
	outTree->Branch("nJets",				&nJets,    				"nJets/I");
	outTree->Branch("nDiJets",				&nDiJets,    			"nDiJets/I");
	outTree->Branch("nBJets",				&nBJets,    			"nBJets/I"); 
	outTree->Branch("nJetsPU",				&nJetsPU,    			"nJetsPU/I"); 
	outTree->Branch("nBJetsPU",				&nBJetsPU,    			"nBJetsPU/I"); 
	outTree->Branch("nLeptons",				&nLeptons, 				"nLeptons/I");  
	outTree->Branch("nLeptons01",			&nLeptons01, 			"nLeptons01/I");
	outTree->Branch("nLeptons04",			&nLeptons04, 			"nLeptons04/I");
	outTree->Branch("matchedBJets",	&matchedBJets,  "matchedBJets/I");
	outTree->Branch("matchedLJets",	&matchedLJets,  "matchedLJets/I");
	outTree->Branch("matchedHiggs",	&matchedHiggs,  "matchedHiggs/I");
	outTree->Branch("bjet1_Pt",		&bjet1_Pt,  	"bjet1_Pt/D");
	outTree->Branch("bjet1_Eta",	&bjet1_Eta,  	"bjet1_Eta/D");
	outTree->Branch("bjet1_Phi",	&bjet1_Phi,  	"bjet1_Phi/D");
	outTree->Branch("bjet1_Mass",	&bjet1_Mass,  	"bjet1_Mass/D");
	outTree->Branch("bjet2_Pt",		&bjet2_Pt,  	"bjet2_Pt/D");
	outTree->Branch("bjet2_Eta",	&bjet2_Eta,  	"bjet2_Eta/D");
	outTree->Branch("bjet2_Phi",	&bjet2_Phi,  	"bjet2_Phi/D");
	outTree->Branch("bjet2_Mass",	&bjet2_Mass,  	"bjet2_Mass/D");
	outTree->Branch("bjet3_Pt",		&bjet3_Pt,  	"bjet3_Pt/D");
	outTree->Branch("bjet3_Eta",	&bjet3_Eta,  	"bjet3_Eta/D");
	outTree->Branch("bjet3_Phi",	&bjet3_Phi,  	"bjet3_Phi/D");
	outTree->Branch("bjet3_Mass",	&bjet3_Mass,  	"bjet3_Mass/D");
	outTree->Branch("bjet4_Pt",		&bjet4_Pt,  	"bjet4_Pt/D");
	outTree->Branch("bjet4_Eta",	&bjet4_Eta,  	"bjet4_Eta/D");
	outTree->Branch("bjet4_Phi",	&bjet4_Phi,  	"bjet4_Phi/D");
	outTree->Branch("bjet4_Mass",	&bjet4_Mass,  	"bjet4_Mass/D");
	outTree->Branch("jet1_Pt",		&jet1_Pt,  		"jet1_Pt/D");
	outTree->Branch("jet1_Eta",		&jet1_Eta,  	"jet1_Eta/D");
	outTree->Branch("jet1_Phi",		&jet1_Phi,  	"jet1_Phi/D");
	outTree->Branch("jet1_Mass",	&jet1_Mass,  	"jet1_Mass/D");
	outTree->Branch("jet2_Pt",		&jet2_Pt,  		"jet2_Pt/D");
	outTree->Branch("jet2_Eta",		&jet2_Eta,  	"jet2_Eta/D");
	outTree->Branch("jet2_Phi",		&jet2_Phi,  	"jet2_Phi/D");
	outTree->Branch("jet2_Mass",	&jet2_Mass,  	"jet2_Mass/D");

  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);
	weight/=Double_t(totalEvents);

    outTree->Fill();

  }
  outFile->Write();
  outFile->Close();
  
	cout << "Finished " << outfilename << " with " << totalEvents << " total events" << endl;

}
