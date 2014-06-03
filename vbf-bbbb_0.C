/*
 * vbf-bbbb_0.C
 * 
 * This macro analyses signal and background samples for vbf-bbbb.
 * It analyses the samples by selecting 4 bjets and 2 light jets.
 * To run use "root vbf-bbbb_0.C"
 * 
 * Code by: Andres Rios
 */

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

// Declare functions

void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void analyzeS();
void analyzeB(TString, Double_t);
void saveResults();
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);
int puJetID(Float_t, Float_t, Float_t);

// Initialize histograms

TH1D *hMassdiJetS = new TH1D("MassdiJetS", "MassdiJetS", 50, 0, 1000);
TH1D *hMass4bJetsS = new TH1D("Mass4bJetsS", "Mass4bJetsS", 100, 0, 1000);
TH1D *hMass2bJetsNoCrossS = new TH1D("Mass2bJetsNoCrossS", "Mass2bJetsNoCrossS", 100, 0, 300);
TH1D *hMass2bJetsCrossS = new TH1D("Mass2bJetsCrossS", "Mass2bJetsCrossS", 100, 0, 300);
TH1D *hInJetsS = new TH1D("InJetsS", "InJetsS", 6, -0.5, 5.5);

TH1D *hMassdiJetB = new TH1D("MassdiJetB", "MassdiJetB", 50, 0, 1000);
TH1D *hMass4bJetsB = new TH1D("Mass4bJetsB", "Mass4bJetsB", 100, 0, 1000);
TH1D *hMass2bJetsNoCrossB = new TH1D("Mass2bJetsNoCrossB", "Mass2bJetsNoCrossB", 100, 0, 300);
TH1D *hMass2bJetsCrossB = new TH1D("Mass2bJetsCrossB", "Mass2bJetsCrossB", 100, 0, 300);
TH1D *hInJetsB = new TH1D("InJetsB", "InJetsB", 6, -0.5, 5.5);

// Initialyze storage variables

Double_t totalSignal = 0, jetsSignal = 0, selJetsSignal = 0, dijetCutSignal = 0, injetCutSignal = 0;
Double_t totalBackground = 0, jetsBackground = 0, selJetsBackground = 0, dijetCutBackground = 0, injetCutBackground = 0;


/*
 * MAIN FUNCTION
 */

void vbf_bbbb_0(){
	
	analyzeS();
	analyzeB("tt-4p-0-600-v1510_14TEV", 530.89358);
	analyzeB("tt-4p-600-1100-v1510_14TEV", 42.55351);
	analyzeB("tt-4p-1100-1700-v1510_14TEV", 4.48209);
	analyzeB("tt-4p-1700-2500-v1510_14TEV", 0.52795);
	analyzeB("tt-4p-2500-100000-v1510_14TEV", 0.05449);
	saveResults();
	
}


/*
 * SIGNAL ANALYSIS
 */

void analyzeS()
{	
	
	const TString inputFilet = "HHToBBBB.txt";
	ifstream ifs(inputFilet);
	assert(ifs.is_open());

	TString filename;
	TChain chain("Delphes");
	
	cout << "Reading signal samples" << endl;

	while(ifs >> filename){
		
		TString filenamef = "root://eoscms.cern.ch//store/user/arapyan/Delphes_phase2/VBFHHTobbbb_TuneZ2_16TeV-madgraph/files/files/" + filename;
		cout << "Reading " << filenamef << endl;
		chain.Add(filenamef);
		
	}
	
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	totalSignal = numberOfEntries;

	// Set up loop variables
	Jet *jet;

	// Set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4, *Jet1, *Jet2;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;
	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vdiJet, v4bJet, v2bJetNoCross, v2bJetCross, vbJet01, vbJet02, vbJet03, vbJet04;
	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");


	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);

	// Reset index variable
	nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;
			
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // Jet loop
			jet = (Jet*) branchJet->At(iJet);
			
			if ((puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0) && (jet->PT > 40) && (fabs(jet->Eta) < 5)){
			
				if(jet->BTag){
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
					
						nbJet3 = nbJet2;
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
				else{
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
		
		
		if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)){
		
			selJetsSignal++;
				
			vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
			vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
					
			vdiJet = vJet1 + vJet2;
				
			hMassdiJetS->Fill(vdiJet.M());
				
			if((vdiJet.M() > 500) && (fabs(Jet1->Eta - Jet2->Eta) > 2.5)){
					
				dijetCutSignal++;
						
				vbJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
				vbJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
				vbJet3.SetPtEtaPhiM(bJet3->PT, bJet3->Eta, bJet3->Phi, bJet3->Mass);
				vbJet4.SetPtEtaPhiM(bJet4->PT, bJet4->Eta, bJet4->Phi, bJet4->Mass);
					
				Int_t inBJets = 0;
				if(Jet1->Eta > Jet2->Eta){
						
					if((bJet1->Eta < Jet1->Eta) && (bJet1->Eta > Jet2->Eta))inBJets++;
					if((bJet2->Eta < Jet1->Eta) && (bJet2->Eta > Jet2->Eta))inBJets++;
					if((bJet3->Eta < Jet1->Eta) && (bJet3->Eta > Jet2->Eta))inBJets++;
					if((bJet4->Eta < Jet1->Eta) && (bJet4->Eta > Jet2->Eta))inBJets++;
						
				}
				else{
						
					if((bJet1->Eta > Jet1->Eta) && (bJet1->Eta < Jet2->Eta))inBJets++;
					if((bJet2->Eta > Jet1->Eta) && (bJet2->Eta < Jet2->Eta))inBJets++;
					if((bJet3->Eta > Jet1->Eta) && (bJet3->Eta < Jet2->Eta))inBJets++;
					if((bJet4->Eta > Jet1->Eta) && (bJet4->Eta < Jet2->Eta))inBJets++;
							
				}
					
				hInJetsS->Fill(inBJets);
					
				if(inBJets==4){
						
					injetCutSignal++;
						
					v4bJet = vbJet1 + vbJet2 + vbJet3 + vbJet4;
					
					hMass4bJetsS->Fill(v4bJet.M());
						
					Jet *bJet01, *bJet02, *bJet03, *bJet04;
						
					bJet04 = bJet1;
					if(bJet2->Eta < bJet04->Eta)bJet04 = bJet2;
					if(bJet3->Eta < bJet04->Eta)bJet04 = bJet3;
					if(bJet4->Eta < bJet04->Eta)bJet04 = bJet4;
							
					bJet01 = bJet04;
					if(bJet1->Eta > bJet01->Eta)bJet01 = bJet1;
					if(bJet2->Eta > bJet01->Eta)bJet01 = bJet2;
					if(bJet3->Eta > bJet01->Eta)bJet01 = bJet3;
					if(bJet4->Eta > bJet01->Eta)bJet01 = bJet4;
							
					bJet02 = bJet04;
					if((bJet1->Eta > bJet02->Eta) && (bJet1 < bJet01))bJet02 = bJet1;
					if((bJet2->Eta > bJet02->Eta) && (bJet2 < bJet01))bJet02 = bJet2;
					if((bJet3->Eta > bJet02->Eta) && (bJet3 < bJet01))bJet02 = bJet3;
					if((bJet4->Eta > bJet02->Eta) && (bJet4 < bJet01))bJet02 = bJet4;
							
					bJet03 = bJet04;
					if((bJet1->Eta > bJet03->Eta) && (bJet1 < bJet02))bJet03 = bJet1;
					if((bJet2->Eta > bJet03->Eta) && (bJet2 < bJet02))bJet03 = bJet2;
					if((bJet3->Eta > bJet03->Eta) && (bJet3 < bJet02))bJet03 = bJet3;
					if((bJet4->Eta > bJet03->Eta) && (bJet4 < bJet02))bJet03 = bJet4;
							
					vbJet01.SetPtEtaPhiM(bJet01->PT, bJet01->Eta, bJet01->Phi, bJet01->Mass);
					vbJet02.SetPtEtaPhiM(bJet02->PT, bJet02->Eta, bJet02->Phi, bJet02->Mass);
					vbJet03.SetPtEtaPhiM(bJet03->PT, bJet03->Eta, bJet03->Phi, bJet03->Mass);
					vbJet04.SetPtEtaPhiM(bJet04->PT, bJet04->Eta, bJet04->Phi, bJet04->Mass);
							
					v2bJetNoCross = vbJet01 + vbJet02;
					hMass2bJetsNoCrossS->Fill(v2bJetNoCross.M());
					v2bJetNoCross = vbJet03 + vbJet04;
					hMass2bJetsNoCrossS->Fill(v2bJetNoCross.M());
							
					v2bJetCross = vbJet01 + vbJet03;
					hMass2bJetsCrossS->Fill(v2bJetCross.M());
					v2bJetCross = vbJet02 + vbJet04;
					hMass2bJetsCrossS->Fill(v2bJetCross.M());
			
				}
				
			}
			
		}
		
	} // End event loop

	
	ifs.close();
}


/*
 * BACKGROUND ANALYSIS
 */

void analyzeB(TString input, Double_t crossSection)
{	
	
	Double_t Background1=0, Background2=0, Background3=0, Background4=0, Background5=0, weight;
	
	const TString inputFilet = input + ".txt";
	ifstream ifs(inputFilet);

	assert(ifs.is_open());

	TString filename;

	TChain chain("Delphes");
	
	cout << "Reading signal samples" << endl;

	while(ifs >> filename){
		TString filenamef = "root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/" + input + "/" + filename;
		cout << "Reading " << filenamef << endl;
		chain.Add(filenamef);
	}
	cout << "Reading " << input << endl;
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Set up loop variables
	Jet *jet;

	// Set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4, *Jet1, *Jet2;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;

	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vdiJet, v4bJet, v2bJetNoCross, v2bJetCross, vbJet01, vbJet02, vbJet03, vbJet04;

	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchEvent = treeReader->UseBranch("Event");

	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);

		// Reset counter variable
		nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;

			
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // Jet loop
			jet = (Jet*) branchJet->At(iJet);
			
			if ((puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0) && (jet->PT > 40) && (fabs(jet->Eta) < 5)){
			
				if(jet->BTag){
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
					
						nbJet3 = nbJet2;
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
				else{
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
		
		
		if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)){
			
			LHEFEvent *event = (LHEFEvent*) branchEvent->At(0);
			weight = event->Weight;
			
			Background2+=weight;
			
			if((Jet1->PT >40) && (Jet2->PT >40) && (bJet1->PT >40) && (bJet2->PT >40) && (bJet3->PT >40) && (bJet4->PT >40) &&
				(fabs(Jet1->Eta)<5) && (fabs(Jet2->Eta)<5) && (fabs(bJet1->Eta)<5) && (fabs(bJet2->Eta)<5) && (fabs(bJet3->Eta)<5) && (fabs(bJet4->Eta)<5)){
				
				Background3+=weight;
				
				vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
				vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
				
				vdiJet = vJet1 + vJet2;
				
				hMassdiJetB->Fill(vdiJet.M());
				
				if((vdiJet.M() > 500) && (fabs(Jet1->Eta - Jet2->Eta) > 2.5)){
					
					Background4+=weight;
					
					Int_t inBJets = 0;
					if(Jet1->Eta > Jet2->Eta){
						
						if((bJet1->Eta < Jet1->Eta) && (bJet1->Eta > Jet2->Eta))inBJets++;
						if((bJet2->Eta < Jet1->Eta) && (bJet2->Eta > Jet2->Eta))inBJets++;
						if((bJet3->Eta < Jet1->Eta) && (bJet3->Eta > Jet2->Eta))inBJets++;
						if((bJet4->Eta < Jet1->Eta) && (bJet4->Eta > Jet2->Eta))inBJets++;
						
					}
					else{
						
						if((bJet1->Eta > Jet1->Eta) && (bJet1->Eta < Jet2->Eta))inBJets++;
						if((bJet2->Eta > Jet1->Eta) && (bJet2->Eta < Jet2->Eta))inBJets++;
						if((bJet3->Eta > Jet1->Eta) && (bJet3->Eta < Jet2->Eta))inBJets++;
						if((bJet4->Eta > Jet1->Eta) && (bJet4->Eta < Jet2->Eta))inBJets++;
						
					}
					
					hInJetsB->Fill(inBJets);
					
					if(inBJets==4){
		
						Background5+=weight;
		
						v4bJet = vbJet1 + vbJet2 + vbJet3 + vbJet4;
					
						hMass4bJetsB->Fill(v4bJet.M());
						
						Jet *bJet01, *bJet02, *bJet03, *bJet04;
							
						bJet04 = bJet1;
						if(bJet2->Eta < bJet04->Eta)bJet04 = bJet2;
						if(bJet3->Eta < bJet04->Eta)bJet04 = bJet3;
						if(bJet4->Eta < bJet04->Eta)bJet04 = bJet4;
						
						bJet01 = bJet04;
						if(bJet1->Eta > bJet01->Eta)bJet01 = bJet1;
						if(bJet2->Eta > bJet01->Eta)bJet01 = bJet2;
						if(bJet3->Eta > bJet01->Eta)bJet01 = bJet3;
						if(bJet4->Eta > bJet01->Eta)bJet01 = bJet4;
						
						bJet02 = bJet04;
						if((bJet1->Eta > bJet02->Eta) && (bJet1 < bJet01))bJet02 = bJet1;
						if((bJet2->Eta > bJet02->Eta) && (bJet2 < bJet01))bJet02 = bJet2;
						if((bJet3->Eta > bJet02->Eta) && (bJet3 < bJet01))bJet02 = bJet3;
						if((bJet4->Eta > bJet02->Eta) && (bJet4 < bJet01))bJet02 = bJet4;
						
						bJet03 = bJet04;
						if((bJet1->Eta > bJet03->Eta) && (bJet1 < bJet02))bJet03 = bJet1;
						if((bJet2->Eta > bJet03->Eta) && (bJet2 < bJet02))bJet03 = bJet2;
						if((bJet3->Eta > bJet03->Eta) && (bJet3 < bJet02))bJet03 = bJet3;
						if((bJet4->Eta > bJet03->Eta) && (bJet4 < bJet02))bJet03 = bJet4;
						
						vbJet01.SetPtEtaPhiM(bJet01->PT, bJet01->Eta, bJet01->Phi, bJet01->Mass);
						vbJet02.SetPtEtaPhiM(bJet02->PT, bJet02->Eta, bJet02->Phi, bJet02->Mass);
						vbJet03.SetPtEtaPhiM(bJet03->PT, bJet03->Eta, bJet03->Phi, bJet03->Mass);
						vbJet04.SetPtEtaPhiM(bJet04->PT, bJet04->Eta, bJet04->Phi, bJet04->Mass);
						
						v2bJetNoCross = vbJet01 + vbJet02;
						hMass2bJetsNoCrossB->Fill(v2bJetNoCross.M());
						v2bJetNoCross = vbJet03 + vbJet04;
						hMass2bJetsNoCrossB->Fill(v2bJetNoCross.M());
						
						v2bJetCross = vbJet01 + vbJet03;
						hMass2bJetsCrossB->Fill(v2bJetCross.M());
						v2bJetCross = vbJet02 + vbJet04;
						hMass2bJetsCrossB->Fill(v2bJetCross.M());
		
					}
					
				}
				
			}
			
		}
		
	} // end event loop
	
	totalBackground += crossSection*3000000*Background1/numberOfEntries;
	jetsBackground += crossSection*3000000*Background2/numberOfEntries;
	selJetsBackground += crossSection*3000000*Background3/numberOfEntries;
	dijetCutBackground += crossSection*3000000*Background4/numberOfEntries;
	injetCutBackground += crossSection*3000000*Background5/numberOfEntries;
	
	ifs.close();
}

/*
 * FUNCTION FOR SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 900, 600);
	
	//gROOT->SetOptStat(kFALSE);

	// Save histograms
	histogram(hMassdiJetS, c1, "Mass diJet system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMassdiJetsS_0.jpg");
	histogram(hMass4bJetsS, c1, "Mass 4 bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass4bJetsS_0.jpg");
	histogram(hMass2bJetsNoCrossS, c1, "Mass 2 not crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsNoCrossS_0.jpg");
	histogram(hMass2bJetsCrossS, c1, "Mass 2 crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsCrossS_0.jpg");
	histogram(hInJetsS, c1, "Number of bjets between light jets", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramInBJetsS_0.jpg");
	
	histogram(hMassdiJetB, c1, "Mass diJet system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMassdiJetsB_0.jpg");
	histogram(hMass4bJetsB, c1, "Mass 4 bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass4bJetsB_0.jpg");
	histogram(hMass2bJetsNoCrossB, c1, "Mass 2 not crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsNoCrossB_0.jpg");
	histogram(hMass2bJetsCrossB, c1, "Mass 2 crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsCrossB_0.jpg");
	histogram(hInJetsB, c1, "Number of bjets between light jets", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramInBJetsB_0.jpg");
	
	// Print event yields
	cout << "Total events: " << 1.8*3000*0.577*0.577 << endl;
	cout << "Events with all jets passing veto: " << 1.8*3000*0.577*0.577*jetsSignal/totalSignal << endl;
	cout << "Events with selected jets: " << 1.8*3000*0.577*0.577*selJetsSignal/totalSignal << endl;
	cout << "Events after dijet cut: " << 1.8*3000*0.577*0.577*dijetCutSignal/totalSignal << endl;
	cout << "Events after injet cut: " << 1.8*3000*0.577*0.577*injetCutSignal/totalSignal << endl;
	
	cout << "\nEvents with all jets passing veto: " << jetsBackground << endl;
	cout << "Events with selected jets: " << selJetsBackground << endl;
	cout << "Events after dijet cut: " << dijetCutBackground << endl;
	cout << "Events after injet cut: " << injetCutBackground << endl;

}

/*
 * FUNCTION FOR SAVING TWO HISTOGRAMS
 */
 
void histogram(TH1D *histoS, TH1D *histoB, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t nS=1, nB=1;
	
	nS/=histoS->Integral();
	nB/=histoB->Integral();
	histoS->Scale(nS);
	histoB->Scale(nB);
	
	Double_t max;
	if((histoS->GetMaximum()) > (histoB->GetMaximum())){max=1.1*(histoS->GetMaximum());}
	else if((histoS->GetMaximum()) <= (histoB->GetMaximum())){max=1.1*(histoB->GetMaximum());}
	
	histoS->SetMaximum(max);
	histoB->SetMaximum(max);
	
	histoS->Draw();
	// add axis labels
	histoS->GetXaxis()->SetTitle(xTitle);
	histoS->GetYaxis()->SetTitle(yTitle);
	histoS->SetTitle(""); // title on top
	
	
	histoB->SetLineColor(kRed);
	histoB->Draw("same");
	

	can->SaveAs(name);
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t norm=1;
	norm/=histo->Integral();
	histo->Scale(norm);
	histo->Draw();
	// add axis labels
	histo->GetXaxis()->SetTitle(xTitle);
	histo->GetYaxis()->SetTitle(yTitle);
	histo->SetTitle(""); // title on top

	can->SaveAs(name);
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

/*
 * FUNCTION FOR JET ID VETO
 */

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
