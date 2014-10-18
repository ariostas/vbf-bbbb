// Include ROOT and C++ libraries
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

// Include Delphes libraries
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

using namespace std;

Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);
int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);


void selectDelphes(const TString inputfile="test.root", const Double_t xsec=0, const Int_t signalFlag=0) {
	
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
		cout << " file broken" << endl;
		return;
	}

	// set up loop variables
	Jet *jet1, *jet2, *jet3, *jet4;
	Photon *photon;
	Electron *electron;
	Muon *muon;
	GenParticle *particle;
	
	Double_t weight=3000000*xsec;

	// set up output variables and file
	UInt_t nEvents, nDiJets, nJets, nBJets, nJetsPU, nBJetsPU, nLeptons, nLeptons01, nLeptons04, matchedBJets, matchedLJets, matchedHiggs;
	
	Double_t bjet1_Pt, bjet1_Eta, bjet1_Phi, bjet1_Mass;
	Double_t bjet2_Pt, bjet2_Eta, bjet2_Phi, bjet2_Mass;
	Double_t bjet3_Pt, bjet3_Eta, bjet3_Phi, bjet3_Mass;
	Double_t bjet4_Pt, bjet4_Eta, bjet4_Phi, bjet4_Mass;
	Double_t jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass;
	Double_t jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass;

	const TString outfile = "/afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp_test/" + inputfile;

	TFile *outFile = new TFile(outfile, "RECREATE");

	// tree to hold the number of events in the file before selection
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
	nEvents=numberOfEntries;
	sampTree->Fill();

	// tree to hold information about selected events
	TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",		&weight,    	"weight/D");  
	outTree->Branch("nJets",		&nJets,    		"nJets/i");
	outTree->Branch("nDiJets",		&nDiJets,    	"nDiJets/i");
	outTree->Branch("nBJets",		&nBJets,    	"nBJets/i"); 
	outTree->Branch("nJetsPU",		&nJetsPU,    	"nJetsPU/i"); 
	outTree->Branch("nBJetsPU",		&nBJetsPU,    	"nBJetsPU/i"); 
	outTree->Branch("nLeptons",		&nLeptons, 		"nLeptons/i");  
	outTree->Branch("nLeptons01",	&nLeptons01, 	"nLeptons01/i");
	outTree->Branch("nLeptons04",	&nLeptons04, 	"nLeptons04/i");
	outTree->Branch("matchedBJets",	&matchedBJets,  "matchedBJets/i");
	outTree->Branch("matchedLJets",	&matchedLJets,  "matchedLJets/i");
	outTree->Branch("matchedHiggs",	&matchedHiggs,  "matchedHiggs/i");
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
	

	for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
		treeReader->ReadEntry(iEntry);
		
		vector<Jet*> selectedJets, orderedJets, selectedLightJets, orderedLightJets;
		TLorentzVector tempv1, tempv2, tempDiJet, tempDiJet1, tempDiJet2;

		// Reset index variables
		nDiJets=nJets=nBJets=nJetsPU=nBJetsPU=nLeptons=nLeptons01=nLeptons04=matchedBJets=matchedLJets=matchedHiggs=0;
		jet1=jet2=jet3=jet4=0;
		selectedJets.clear(); orderedJets.clear(); selectedLightJets.clear(); orderedLightJets.clear();
		
		weight=3000000*xsec;
		if(signalFlag == 0) weight *= ((LHEFEvent*) branchEvent->At(0))->Weight;
		
		for (Int_t i=0; i<branchJet->GetEntries(); i++) {
			
			Double_t highestPt=-1, highestPtLight=-1;
			Jet *highestPtJet=0, *highestPtLightJet=0;
			
			for (Int_t j=0; j<branchJet->GetEntries(); j++) {
				jet1 = (Jet*) branchJet->At(j);
				
				if((jet1->PT <= 30) || (fabs(jet1->Eta) >= 4)) continue;
			
				bool count = true;
					
				for (Int_t iP=0; iP<branchElectron->GetEntries(); iP++) {
					electron = (Electron*) branchElectron->At(iP);
					if((electron->PT > 25) && (electron->IsolationVar < 0.4) && (deltaR(electron->Eta, jet1->Eta, electron->Phi, jet1->Phi) < 0.4)) count = false;
				}
				for (Int_t iP=0; iP<branchMuon->GetEntries(); iP++) {
					muon = (Muon*) branchMuon->At(iP);
					if((muon->PT > 25) && (muon->IsolationVar < 0.4) && (deltaR(muon->Eta, jet1->Eta, muon->Phi, jet1->Phi) < 0.4)) count = false;
				}
				for (Int_t iP=0; iP<branchPhoton->GetEntries(); iP++) {
					photon = (Photon*) branchPhoton->At(iP);
					if((photon->PT > 25) && (photon->IsolationVar < 0.4) && (deltaR(photon->Eta, jet1->Eta, photon->Phi, jet1->Phi) < 0.4)) count = false;
				}
			
				//Ignore jets near isolated objects
				if(!count) continue;
				
				if ((i==0) && (fabs(jet1->Eta)<4.0) && (jet1->PT > 30)){
				
					nJets++;
					if (jet1->BTag > 1) nBJets++;
					
					if(puJetID(jet1->Eta, jet1->MeanSqDeltaR, jet1->BetaStar) == 0){
						
						nJetsPU++;
						if (jet1->BTag > 1) nBJetsPU++;
					}
				}
				
				if(puJetID(jet1->Eta, jet1->MeanSqDeltaR, jet1->BetaStar) != 0) continue;
				
				if((jet1->BTag >= 2) && (fabs(jet1->Eta) < 2.5)){
				
					bool isDifferent = true;
					for(UInt_t x=0; x<orderedJets.size(); x++){
						if(jet1 == orderedJets.at(x)) isDifferent = false;
					}
					if(!isDifferent || jet1->PT <= highestPt) continue;
				
					highestPt = jet1->PT;
					highestPtJet = jet1;
					
				}
				
				else if(jet1->BTag < 2){
				
					bool isDifferent = true;
					for(UInt_t x=0; x<orderedLightJets.size(); x++){
						if(jet1 == orderedLightJets.at(x)) isDifferent = false;
					}
					if(!isDifferent || jet1->PT <= highestPtLight) continue;
				
					highestPtLight = jet1->PT;
					highestPtLightJet = jet1;
					
				}
				
			}
			
			if(highestPt != -1 && highestPtJet) orderedJets.push_back(highestPtJet);
			
			if(highestPtLight != -1 && highestPtLightJet) orderedLightJets.push_back(highestPtLightJet);	
			
		}
		
		if((orderedJets.size() < 4) || (orderedLightJets.size() < 2)) continue;
		
		// Select bjets
		for (UInt_t i=0; i<orderedJets.size(); i++) { // bjet loop
			jet1 = orderedJets.at(i);
				
			tempv1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
			
			for (UInt_t j=0; j<orderedJets.size(); j++) {
				jet2 = orderedJets.at(j);
				
				if(i==j) continue;
					
				bool isDifferent = true;
				for(UInt_t x=0; x<selectedJets.size(); x++){
					if(jet1 == selectedJets.at(x) || jet2 == selectedJets.at(x)) isDifferent = false;
				}
				if(!isDifferent) continue;
				
				Double_t dRtemp = deltaR(jet1->Eta, jet2->Eta, jet1->Phi, jet2->Phi);
					
				tempv2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
				tempDiJet = tempv1 + tempv2;
				
				if((tempDiJet.Pt() <= 100) || (tempDiJet.M() <= 80) || (tempDiJet.M() >= 140) || (dRtemp <= 0.5) || (dRtemp >= 2)) continue;
					
				selectedJets.push_back(jet1);
				selectedJets.push_back(jet2);
				
			}
			
		}// End bjet loop
		
		orderedJets.clear();
		
		if(selectedJets.size()%2 != 0){cout << "\n\n\n\n\nError. Check the code\n\n\n\n\n" << endl; return;}
		
		nDiJets = selectedJets.size()/2;
		
		if(nDiJets < 2) continue;
		
		jet1 = selectedJets.at(0);
		jet2 = selectedJets.at(1);
		jet3 = selectedJets.at(2);
		jet4 = selectedJets.at(3);
		
		tempv1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
		tempv2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
		
		tempDiJet1 = tempv1 + tempv2;
		
		tempv1.SetPtEtaPhiM(jet3->PT, jet3->Eta, jet3->Phi, jet3->Mass);
		tempv2.SetPtEtaPhiM(jet4->PT, jet4->Eta, jet4->Phi, jet4->Mass);
		
		tempDiJet2 = tempv1 + tempv2;
		
		// Select light jets
		for (UInt_t i=0; i<orderedLightJets.size(); i++) { // jet loop
			jet1 = orderedLightJets.at(i);
				
			tempv1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
			
			for (UInt_t j=0; j<orderedLightJets.size(); j++) {
				jet2 = orderedLightJets.at(j);
				
				if(i==j) continue;
					
				bool isDifferent = true;
				for(UInt_t x=0; x<selectedLightJets.size(); x++){
					if(jet1 == selectedLightJets.at(x) || jet2 == selectedLightJets.at(x)) isDifferent = false;
				}
				if(!isDifferent) continue;
				
				Double_t dRtemp = deltaR(jet1->Eta, jet2->Eta, jet1->Phi, jet2->Phi);
					
				tempv2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
				tempDiJet = tempv1 + tempv2;
				
				if((tempDiJet.Pt() <= 50) || (tempDiJet.M() <= 250) || (fabs(jet1->Eta - jet2->Eta) <= 3.5)) continue;
				
				Int_t nInnerDiJets = 0;
				if(jet1->Eta > jet2->Eta){
			
					if((tempDiJet1.Eta() < jet1->Eta) && (tempDiJet1.Eta() > jet2->Eta))nInnerDiJets++;
					if((tempDiJet2.Eta() < jet1->Eta) && (tempDiJet2.Eta() > jet2->Eta))nInnerDiJets++;
			
				}
				else{
			
					if((tempDiJet1.Eta() > jet1->Eta) && (tempDiJet1.Eta() < jet2->Eta))nInnerDiJets++;
					if((tempDiJet2.Eta() > jet1->Eta) && (tempDiJet2.Eta() < jet2->Eta))nInnerDiJets++;
				
				}
				
				if(nInnerDiJets != 2) continue;
					
				selectedLightJets.push_back(jet1);
				selectedLightJets.push_back(jet2);
				
			}
			
		}// End jet loop
		
		orderedLightJets.clear();
		
		if(selectedLightJets.size()%2 != 0){cout << "\n\n\n\n\nError. Check the code\n\n\n\n\n" << endl; return;}
		
		if(selectedLightJets.size() < 2) continue;
		
		bool bj1isMatched=false, bj2isMatched=false, bj3isMatched=false, bj4isMatched=false, j1isMatched=false, j2isMatched=false, h1isMatched=false, h2isMatched=false;
		
		for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { 
			particle = (GenParticle*) branchParticle->At(iParticle);
			
			if(abs(particle->PID) == 5){
				
				if(deltaR(particle->Eta, selectedJets.at(0)->Eta, particle->Phi, selectedJets.at(0)->Phi) < 0.4 && !bj1isMatched){ matchedBJets++; bj1isMatched=true;}
				else if(deltaR(particle->Eta, selectedJets.at(1)->Eta, particle->Phi, selectedJets.at(1)->Phi) < 0.4 && !bj2isMatched){ matchedBJets++; bj2isMatched=true;}
				else if(deltaR(particle->Eta, selectedJets.at(2)->Eta, particle->Phi, selectedJets.at(2)->Phi) < 0.4 && !bj3isMatched){ matchedBJets++; bj3isMatched=true;}
				else if(deltaR(particle->Eta, selectedJets.at(3)->Eta, particle->Phi, selectedJets.at(3)->Phi) < 0.4 && !bj4isMatched){ matchedBJets++; bj4isMatched=true;}
				
			}
			
			else if(abs(particle->PID) >= 1 && abs(particle->PID) <= 4){
				
				if(deltaR(particle->Eta, selectedLightJets.at(0)->Eta, particle->Phi, selectedLightJets.at(0)->Phi) < 0.4 && !j1isMatched){ matchedLJets++; j1isMatched=true;}
				else if(deltaR(particle->Eta, selectedLightJets.at(1)->Eta, particle->Phi, selectedLightJets.at(1)->Phi) < 0.4 && !j2isMatched){ matchedLJets++; j2isMatched=true;}
				
			}
			
			if(abs(particle->PID) == 25){
				
				if(deltaR(particle->Eta, tempDiJet1.Eta(), particle->Phi, tempDiJet1.Phi()) < 0.4 && !h1isMatched){ matchedHiggs++; h1isMatched=true;}
				else if(deltaR(particle->Eta, tempDiJet2.Eta(), particle->Phi, tempDiJet2.Phi()) < 0.4 && !h2isMatched){ matchedHiggs++; h2isMatched=true;}
				
			}
		}
			
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
		
		bjet1_Pt = 		selectedJets.at(0)->PT;
		bjet1_Eta = 	selectedJets.at(0)->Eta;
		bjet1_Phi = 	selectedJets.at(0)->Phi;
		bjet1_Mass = 	selectedJets.at(0)->Mass;
		
		bjet2_Pt = 		selectedJets.at(1)->PT;
		bjet2_Eta = 	selectedJets.at(1)->Eta;
		bjet2_Phi = 	selectedJets.at(1)->Phi;
		bjet2_Mass = 	selectedJets.at(1)->Mass;
		
		bjet3_Pt = 		selectedJets.at(2)->PT;
		bjet3_Eta = 	selectedJets.at(2)->Eta;
		bjet3_Phi = 	selectedJets.at(2)->Phi;
		bjet3_Mass = 	selectedJets.at(2)->Mass;
		
		bjet4_Pt = 		selectedJets.at(3)->PT;
		bjet4_Eta = 	selectedJets.at(3)->Eta;
		bjet4_Phi = 	selectedJets.at(3)->Phi;
		bjet4_Mass = 	selectedJets.at(3)->Mass;
		
		jet1_Pt = 	selectedLightJets.at(0)->PT;
		jet1_Eta = 	selectedLightJets.at(0)->Eta;
		jet1_Phi = 	selectedLightJets.at(0)->Phi;
		jet1_Mass = selectedLightJets.at(0)->Mass;
		
		jet2_Pt = 	selectedLightJets.at(1)->PT;
		jet2_Eta = 	selectedLightJets.at(1)->Eta;
		jet2_Phi = 	selectedLightJets.at(1)->Phi;
		jet2_Mass = selectedLightJets.at(1)->Mass;
		
		outTree->Fill();
		

	} // end event loop

	// save file
	outFile->Write();
	// close file
	outFile->Close();

	cout << " selection done." << endl;

	

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

/*
 * FUNCTION TO CALCULATE DELTA R
 */

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

	const Float_t pi = 3.14159265358979;

	Float_t etaDiff = (eta1-eta2);
	Float_t phiDiff = fabs(phi1-phi2);
	while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

	Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

	return TMath::Sqrt(deltaRSquared);

}
