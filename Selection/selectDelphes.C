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

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void selectDelphes(const TString inputfile="test.root", const Double_t xsec=0, const Int_t signalFlag=0){
	
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
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	
	TString tempString = "Reading " + inputfile + "...";
	tempString.Resize(75);
	cout << tempString;
	
	if (!(branchJet)) {
		cout << "  file broken" << endl;
		return;
	}

	// set up loop variables
	Jet *jet;
	Electron *electron;
	Photon *photon;
	Muon *muon;
	GenParticle *particle;

	// set up storage variables
	Jet *bJet1=0, *bJet2=0, *bJet3=0, *bJet4=0, *Jet1=0, *Jet2=0;
	Double_t weight=3000000*xsec;
	UInt_t nEvents=0;
	Int_t nLeptons=0, nLeptons01=0, nLeptons04=0, nJets=0, nJetsNoPU=0;
	Int_t matchedBJets=0, matchedLJets=0;
	
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
	outTree->Branch("nJetsNoPU",	&nJetsNoPU,  	"nJetsNoPU/I");
	outTree->Branch("matchedBJets",	&matchedBJets,  "matchedBJets/I");
	outTree->Branch("matchedLJets",	&matchedLJets,  "matchedLJets/I");


	for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
		treeReader->ReadEntry(iEntry);

		// Reset variables
		nLeptons=nLeptons01=nLeptons04=nJets=nJetsNoPU=matchedBJets=matchedLJets=0;
		bJet1=bJet2=bJet3=bJet4=Jet1=Jet2=0;
		
		weight=3000000*xsec;
		if(signalFlag == 0) weight *= ((LHEFEvent*) branchEvent->At(0))->Weight;
		
		// Select bjets and light jets
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // Jet loop
			jet = (Jet*) branchJet->At(iJet);
			
			if(jet->PT <= 30 || jet->Eta >= 4) continue;
			
			bool count = true;
					
			for (Int_t iP=0; iP<branchElectron->GetEntries(); iP++) {
				electron = (Electron*) branchElectron->At(iP);
				if(electron->PT > 25 && electron->IsolationVar < 0.4 && deltaR(electron->Eta, jet->Eta, electron->Phi, jet->Phi) < 0.4) count = false;
			}
			for (Int_t iP=0; iP<branchMuon->GetEntries(); iP++) {
				muon = (Muon*) branchMuon->At(iP);
				if(muon->PT > 25 && muon->IsolationVar < 0.4 && deltaR(muon->Eta, jet->Eta, muon->Phi, jet->Phi) < 0.4) count = false;
			}
			for (Int_t iP=0; iP<branchPhoton->GetEntries(); iP++) {
				photon = (Photon*) branchPhoton->At(iP);
				if(photon->PT > 25 && photon->IsolationVar < 0.4 && deltaR(photon->Eta, jet->Eta, photon->Phi, jet->Phi) < 0.4) count = false;
			}
			
			//Ignore jets near isolated objects
			if(!count) continue;
			
			nJets++;

			if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0){
				
				nJetsNoPU++;

				if((jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7) && (fabs(jet->Eta) < 2.5)){
					
					if(!bJet1){
						
						bJet1 = (Jet*) branchJet->At(iJet);
						
					}
					else if(jet->PT > bJet1->PT){

						bJet4 = bJet3;
						bJet3 = bJet2;
						bJet2 = bJet1;
						bJet1 = (Jet*) branchJet->At(iJet);

					}
					else if(!bJet2){

						bJet2 = (Jet*) branchJet->At(iJet);

					}
					else if(jet->PT > bJet2->PT){

						bJet4 = bJet3;
						bJet3 = bJet2;
						bJet2 = (Jet*) branchJet->At(iJet);

					}
					else if(!bJet3){

						bJet3 = (Jet*) branchJet->At(iJet);

					}
					else if(jet->PT > bJet3->PT){

						bJet4 = bJet3;
						bJet3 = (Jet*) branchJet->At(iJet);

					}
					else if(!bJet4){

						bJet4 = (Jet*) branchJet->At(iJet);

					}
					else if(jet->PT > bJet4->PT){

						bJet4 = (Jet*) branchJet->At(iJet);

					}
				}
				else if(!(jet->BTag == 2 || jet->BTag == 3 || jet->BTag == 6 || jet->BTag == 7)){
					if(!Jet1){

						Jet1 = (Jet*) branchJet->At(iJet);

					}
					else if((jet->PT > Jet1->PT)){

						Jet2 = Jet1;
						Jet1 = (Jet*) branchJet->At(iJet);

					}
					else if(!Jet2){

						Jet2 = (Jet*) branchJet->At(iJet);

					}
					else if(jet->PT > Jet2->PT){

						Jet2 = (Jet*) branchJet->At(iJet);

					}

				}

			}

		} // End jet loop
		
		// Check if there are four bjets and two light jets
		if((bJet1) && (bJet2) && (bJet3) && (bJet4) && (Jet1) && (Jet2)){
			
			// Match with generator level objects
			bool bj1isMatched=false, bj2isMatched=false, bj3isMatched=false, bj4isMatched=false, j1isMatched=false, j2isMatched=false;
			
			for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { 
				particle = (GenParticle*) branchParticle->At(iParticle);
				
				if(abs(particle->PID) == 5){
				
				if(deltaR(particle->Eta, bJet1->Eta, particle->Phi, bJet1->Phi) < 0.4 && !bj1isMatched){ matchedBJets++; bj1isMatched=true;}
				else if(deltaR(particle->Eta, bJet2->Eta, particle->Phi, bJet2->Phi) < 0.4 && !bj2isMatched){ matchedBJets++; bj2isMatched=true;}
				else if(deltaR(particle->Eta, bJet3->Eta, particle->Phi, bJet3->Phi) < 0.4 && !bj3isMatched){ matchedBJets++; bj3isMatched=true;}
				else if(deltaR(particle->Eta, bJet4->Eta, particle->Phi, bJet4->Phi) < 0.4 && !bj4isMatched){ matchedBJets++; bj4isMatched=true;}
				
				}
			
				else if(abs(particle->PID) >= 1 && abs(particle->PID) <= 4){
				
					if(deltaR(particle->Eta, Jet1->Eta, particle->Phi, Jet1->Phi) < 0.4 && !j1isMatched){ matchedLJets++; j1isMatched=true;}
					else if(deltaR(particle->Eta, Jet2->Eta, particle->Phi, Jet2->Phi) < 0.4 && !j2isMatched){ matchedLJets++; j2isMatched=true;}
				
				}
			}
			
			// Count leptons
			nLeptons = branchElectron->GetEntries() + branchMuon->GetEntries();
			
			for (Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++) { 
				electron = (Electron*) branchElectron->At(iElectron);
				if(electron->IsolationVar < 0.1 && electron->PT > 25) nLeptons01++;
				if(electron->IsolationVar < 0.4 && electron->PT > 25) nLeptons04++;
			}
			
			for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { 
				muon = (Muon*) branchMuon->At(iMuon);
				if(muon->IsolationVar < 0.1 && muon->PT > 25) nLeptons01++;
				if(muon->IsolationVar < 0.4 && muon->PT > 25) nLeptons04++;
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


Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar) {
  
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
