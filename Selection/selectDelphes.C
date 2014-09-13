//-------------------------------------------------------------------
// Select vbf-bbbb events
//-------------------------------------------------------------------

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

#include "../Delphes/modules/Delphes.h"                   // delphes
#include "../Delphes/external/ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "../Delphes/classes/DelphesClasses.h"            // delphes

#endif

using namespace std;

int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void selectDelphes(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToGGBB_14TeV/HHToGGBB_14TeV_0.root",
const Double_t xsec=1341.36923, Int_t signalFlag=1) {
	
	TString inputFile = inputfile;
	if(signalFlag == 1) inputFile.Remove(0,15);

	// read input input file
	TChain chain("Delphes");
	chain.Add("root://eoscms.cern.ch/" + (signalFlag == 1 ? "/store/group/upgrade/dihiggs_signal_4b/VBFHHTobbbb_TuneZ2_14TeV-madgraph/files/" /*"/store/group/upgrade/delphes/dihiggs_signal_bbbb/gFHHTobbbb_TuneZ2_14TeV_madgraph/files/"*/  + inputFile: "/store/group/upgrade/delphes/PhaseII_140PU_ProdJul28/" + inputFile));
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// set up branches to read in from file
	TClonesArray *branchEvent;
	if(signalFlag == 0) branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	
	TString tempString = "Reading " + inputfile + "...";
	tempString.Resize(75);
	cout << tempString;
	
	if (!(branchJet)) {
		cout << " file broken" << endl;
		return;
	}

	// set up loop variables
	Jet *jet;
	Electron *electron;
	Photon *photon;
	Muon *muon;

	// set up storage variables
	Jet *bJet1=0, *bJet2=0, *bJet3=0, *bJet4=0, *Jet1=0, *Jet2=0;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;
	Double_t weight=3000000*xsec;
	UInt_t nEvents=0;
	Int_t nLeptons=0, nLeptons01=0, nLeptons04=0, nJets=0, nJetsHighPt=0, nCentralLightJets=0, nCentralLightJetsPU=0;
	
	Double_t bjet1_Pt, bjet1_Eta, bjet1_Phi, bjet1_Mass;
	Double_t bjet2_Pt, bjet2_Eta, bjet2_Phi, bjet2_Mass;
	Double_t bjet3_Pt, bjet3_Eta, bjet3_Phi, bjet3_Mass;
	Double_t bjet4_Pt, bjet4_Eta, bjet4_Phi, bjet4_Mass;
	Double_t jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass;
	Double_t jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass;

	const TString outfile = "/afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/" + inputfile;

	TFile *outFile = new TFile(outfile, "RECREATE");

	// tree to hold the number of events in the file before selection
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
	nEvents=numberOfEntries;
	sampTree->Fill();

	// tree to hold information about selected events
	TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",		&weight,    	"weight/D");
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
	outTree->Branch("nLeptons",		&nLeptons,    	"nLeptons/I");
	outTree->Branch("nLeptons01",	&nLeptons01,   	"nLeptons01/I");
	outTree->Branch("nLeptons04",	&nLeptons04,   	"nLeptons04/I");
	outTree->Branch("nJets",		&nJets,    		"nJets/I");
	outTree->Branch("nJetsHighPt",	&nJetsHighPt,  	"nJetsHighPt/I");
	outTree->Branch("nCentralLightJets",	&nCentralLightJets,  	"nCentralLightJets/I");
	outTree->Branch("nCentralLightJetsPU",	&nCentralLightJetsPU,  	"nCentralLightJetsPU/I");


	for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
		treeReader->ReadEntry(iEntry);

		// Reset index and storage variables
		nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;
		nLeptons=nLeptons01=nLeptons04=nJets=nJetsHighPt=nCentralLightJets=nCentralLightJetsPU=0;
		
		weight=3000000*xsec;
		if(signalFlag == 0) weight *= ((LHEFEvent*) branchEvent->At(0))->Weight;
		
		// Select bjets and light jets
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // Jet loop
			jet = (Jet*) branchJet->At(iJet);
			
			if(jet->PT <= 30) continue;
			
			bool count = true;
					
			for (Int_t iP=0; iP<branchElectron->GetEntries(); iP++) {
				electron = (Electron*) branchElectron->At(iP);
				if(electron->PT > 10 && electron->IsolationVar < 0.4 && deltaR(electron->Eta, jet->Eta, electron->Phi, jet->Phi) < 0.4) count = false;
			}
			for (Int_t iP=0; iP<branchMuon->GetEntries(); iP++) {
				muon = (Muon*) branchMuon->At(iP);
				if(muon->PT > 10 && muon->IsolationVar < 0.4 && deltaR(muon->Eta, jet->Eta, muon->Phi, jet->Phi) < 0.4) count = false;
			}
			for (Int_t iP=0; iP<branchPhoton->GetEntries(); iP++) {
				photon = (Photon*) branchPhoton->At(iP);
				if(photon->PT > 10 && photon->IsolationVar < 0.4 && deltaR(photon->Eta, jet->Eta, photon->Phi, jet->Phi) < 0.4) count = false;
			}
			
			//Ignore jets near isolated objects
			if(!count) continue;
			
			if((jet->PT > 30) && (fabs(jet->Eta) < 1.0) && !(jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7)) nCentralLightJets++;

			if ((puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0) && (jet->PT > 30) && (fabs(jet->Eta) < 4.0)){
				
				nJets++;
				if((fabs(jet->Eta) < 1.0) && !(jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7)) nCentralLightJetsPU++;
				if(jet->PT > 40) nJetsHighPt++;

				if((jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7) && (fabs(jet->Eta) < 2.5)){
					
					if(nbJet1 == -1){
						
						nbJet1 = iJet;
						bJet1 = (Jet*) branchJet->At(nbJet1);
						
					}
					else if(jet->PT > bJet1->PT){

						nbJet4 = nbJet3;
						bJet4 = (Jet*) branchJet->At(nbJet4);

						nbJet3 = nbJet2;
						bJet3 = (Jet*) branchJet->At(nbJet3);

						nbJet2 = nbJet1;
						bJet2 = (Jet*) branchJet->At(nbJet2);

						nbJet1 = iJet;
						bJet1 = (Jet*) branchJet->At(nbJet1);

					}
					else if(nbJet2 == -1){

						nbJet2 = iJet;
						bJet2 = (Jet*) branchJet->At(nbJet2);

					}
					else if(jet->PT > bJet2->PT){

						nbJet4 = nbJet3;
						bJet4 = (Jet*) branchJet->At(nbJet4);

						nbJet3 = nbJet2;
						bJet3 = (Jet*) branchJet->At(nbJet3);

						nbJet2 = iJet;
						bJet2 = (Jet*) branchJet->At(nbJet2);

					}
					else if(nbJet3 == -1){

					nbJet3 = iJet;
					bJet3 = (Jet*) branchJet->At(nbJet3);

					}
					else if(jet->PT > bJet3->PT){

						nbJet4 = nbJet3;
						bJet4 = (Jet*) branchJet->At(nbJet4);

						nbJet3 = iJet;
						bJet3 = (Jet*) branchJet->At(nbJet3);

					}
					else if(nbJet4 == -1){

						nbJet4 = iJet;
						bJet4 = (Jet*) branchJet->At(nbJet4);

					}
					else if(jet->PT > bJet4->PT){

						nbJet4 = iJet;
						bJet4 = (Jet*) branchJet->At(nbJet4);

					}
				}
				else if(!(jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7)){
					if(nJet1 == -1){

						nJet1 = iJet;
						Jet1 = (Jet*) branchJet->At(nJet1);

					}
					else if((jet->PT > Jet1->PT)){

						nJet2 = nJet1;
						Jet2 = (Jet*) branchJet->At(nJet2);

						nJet1 = iJet;
						Jet1 = (Jet*) branchJet->At(nJet1);

					}
					else if(nJet2 == -1){

						nJet2 = iJet;
						Jet2 = (Jet*) branchJet->At(nJet2);

					}
					else if(jet->PT > Jet2->PT){

						nJet2 = iJet;
						Jet2 = (Jet*) branchJet->At(nJet2);

					}

				}

			}

		} // End jet loop
		
		// Check if there are four bjets and two light jets
		if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)){
			
			nLeptons = branchElectron->GetEntries() + branchMuon->GetEntries();
			
			for (Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++) { 
				electron = (Electron*) branchElectron->At(iElectron);
				if(electron->IsolationVar < 0.1 && electron->PT > 10) nLeptons01++;
				if(electron->IsolationVar < 0.4 && electron->PT > 10) nLeptons04++;
			}
			
			for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { 
				muon = (Muon*) branchMuon->At(iMuon);
				if(muon->IsolationVar < 0.1 && muon->PT > 10) nLeptons01++;
				if(muon->IsolationVar < 0.4 && muon->PT > 10) nLeptons04++;
			}
			
			bjet1_Pt = bJet1->PT;
			bjet1_Eta = bJet1->Eta;
			bjet1_Phi = bJet1->Phi;
			bjet1_Mass = bJet1->Mass;
			
			bjet2_Pt = bJet2->PT;
			bjet2_Eta = bJet2->Eta;
			bjet2_Phi = bJet2->Phi;
			bjet2_Mass = bJet2->Mass;
			
			bjet3_Pt = bJet3->PT;
			bjet3_Eta = bJet3->Eta;
			bjet3_Phi = bJet3->Phi;
			bjet3_Mass = bJet3->Mass;
			
			bjet4_Pt = bJet4->PT;
			bjet4_Eta = bJet4->Eta;
			bjet4_Phi = bJet4->Phi;
			bjet4_Mass = bJet4->Mass;
			
			jet1_Pt = Jet1->PT;
			jet1_Eta = Jet1->Eta;
			jet1_Phi = Jet1->Phi;
			jet1_Mass = Jet1->Mass;
			
			jet2_Pt = Jet2->PT;
			jet2_Eta = Jet2->Eta;
			jet2_Phi = Jet2->Phi;
			jet2_Mass = Jet2->Mass;
			
			outTree->Fill();
			
		}

	} // end event loop

	// save file
	outFile->Write();
	// close file
	outFile->Close();

	cout << "  selection done." << endl;

	

}


int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar) {
  
	Float_t MeanSqDeltaRMaxBarrel=0.07;
	Float_t BetaMinBarrel=0.87;
	Float_t MeanSqDeltaRMaxEndcap=0.07;
	Float_t BetaMinEndcap=0.85;

	//cout << eta << ", " << meanSqDeltaR << ", " << betastar << ": ";

	if (fabs(eta)<1.5) {
		if ((meanSqDeltaR<MeanSqDeltaRMaxBarrel)&&(betastar<BetaMinBarrel)) {
			//cout << "barrel 0" << endl;
			return 0;
		}
		else {
			//cout << "barrel 1" << endl;
			return 1;
		}
	}
	else if (fabs(eta)<4.0) {
		if ((meanSqDeltaR<MeanSqDeltaRMaxEndcap)&&(betastar<BetaMinEndcap)) {
			//cout << "endcap 0" << endl;
			return 0;
		}
		else {
			//cout << "endcap 1" << endl;
			return 1;
		}
	}
	//cout << "forward 1" << endl;
	return 1;

}

Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
