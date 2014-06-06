/*
 * vbf-bbbb_3.C
 * 
 * This macro analyses signal and background samples for vbf-bbbb.
 * It analyses the samples by selecting 6 jets, and taking the leading jets to be the light jets.
 * To run use "root vbf-bbbb_3.C"
 * 
 * The selection of events consist in requiring four bjets and two light
 * jets that pass the jet id veto, have pt > 40 and |eta| < 5
 * 
 * The dijet cut consist in requiring the light jets to yield a reconstructed
 * mass > 500 GeV and to have delta eta > 2.5
 * 
 * The injet cut requires all bjets to be between the light jets.
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
#include <TStyle.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

// Include Delphes libraries
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

// Declare constant variables
const bool Signal=true, Background=false;

// Declare functions
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void analyze(TString, Double_t, bool);
void saveResults();
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
Double_t selectionSignal = 0, dijetCutSignal = 0, injetCutSignal = 0;
Double_t selectionBackground = 0, dijetCutBackground = 0, injetCutBackground = 0;



/*
 * MAIN FUNCTION
 */

void vbf_bbbb_3(){
	
	// Analyze signal
	analyze("HHToBBBB_14TeV", 0.599, Signal);
	
	// Analyse ttbar background
	analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background);
	analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background);
	analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background);
	analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background);
	analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background);
	
	// Save results
	saveResults();
	
}


/*
 * ANALYSIS
 */

void analyze(TString inputfile, Double_t crossSection, bool SorB)
{	
	// The SorB variable is true if it's signal and false if it's a background
	
	const TString inputFile = inputfile + ".txt";
	ifstream ifs(inputFile);

	assert(ifs.is_open());

	TString filename;

	TChain chain("Delphes");
	
	cout << "Reading " << inputfile << endl;

	while(ifs >> filename){
		
		TString filenamef;
		
		if(SorB) filenamef = "root://eoscms.cern.ch//store/user/arapyan/Delphes_phase2/VBFHHTobbbb_TuneZ2_16TeV-madgraph/files/files/" + filename;
		else filenamef = "root://eoscms.cern.ch//eos/cms/store/group/upgrade/delphes/test4/" + filename;
		
		cout << "Reading " << filenamef << endl;
		chain.Add(filenamef);
		
	}
	
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Set up loop variables
	Jet *jet;

	// Set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4, *Jet1, *Jet2, *Jet3, *Jet4, *Jet5, *Jet6;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2, nJet3, nJet4, nJet5, nJet6;
	LHEFEvent *event;
	Double_t weight;
	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vdiJet, v4bJet, v2bJetNoCross, v2bJetCross, vbJet01, vbJet02, vbJet03, vbJet04;

	// Set up temporal variables
	Double_t tempSelection=0, tempDijetCut=0, tempInjetCut=0;
	
	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchEvent;
	(!SorB ? branchEvent = treeReader->UseBranch("Event") : branchEvent = 0);
	

	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // Event loop
		treeReader->ReadEntry(iEntry);

		// Reset index variables
		nJet1=nJet2=nJet3=nJet4=nJet5=nJet6=-1;

			
		for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // Jet loop
			jet = (Jet*) branchJet->At(iJet);
			
			if ((puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar) == 0) && (jet->PT > 40) && (fabs(jet->Eta) < 5)){
				
					if(nJet1 == -1){
					
					nJet1 = iJet;
					Jet1 = (Jet*) branchJet->At(nJet1);
					
				}
				else if((jet->PT > Jet1->PT)){
					
					nJet6 = nJet5;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
					nJet5 = nJet4;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
					nJet4 = nJet3;
					Jet4 = (Jet*) branchJet->At(nJet4);
					
					nJet3 = nJet2;
					Jet3 = (Jet*) branchJet->At(nJet3);
					
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
					
					nJet6 = nJet5;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
					nJet5 = nJet4;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
					nJet4 = nJet3;
					Jet4 = (Jet*) branchJet->At(nJet4);
					
					nJet3 = nJet2;
					Jet3 = (Jet*) branchJet->At(nJet3);
					
					nJet2 = iJet;
					Jet2 = (Jet*) branchJet->At(nJet2);
					
				}
				else if(nJet3==-1){
					
					nJet3 = iJet;
					Jet3 = (Jet*) branchJet->At(nJet3);
					
				}
				else if(jet->PT > Jet3->PT){
					
					nJet6 = nJet5;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
					nJet5 = nJet4;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
					nJet4 = nJet3;
					Jet4 = (Jet*) branchJet->At(nJet4);
					
					nJet3 = iJet;
					Jet3 = (Jet*) branchJet->At(nJet3);
					
				}
				else if(nJet4==-1){
					
					nJet4 = iJet;
					Jet4 = (Jet*) branchJet->At(nJet4);
					
				}
				else if(jet->PT > Jet4->PT){
					
					nJet6 = nJet5;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
					nJet5 = nJet4;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
					nJet4 = iJet;
					Jet4 = (Jet*) branchJet->At(nJet4);
					
				}
				else if(nJet5==-1){
					
					nJet5 = iJet;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
				}
				else if(jet->PT > Jet5->PT){
					
					nJet6 = nJet5;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
					nJet5 = iJet;
					Jet5 = (Jet*) branchJet->At(nJet5);
					
				}
				else if(nJet6==-1){
					
					nJet6 = iJet;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
				}
				else if(jet->PT > Jet6->PT){
					
					nJet6 = iJet;
					Jet6 = (Jet*) branchJet->At(nJet6);
					
				}
			
			}
			
		}// End jet loop
		
		// Check for six selected jets
		if((nJet1!=-1) && (nJet2!=-1) && (nJet3!=-1) && (nJet4!=-1) && (nJet5!=-1) && (nJet6!=-1)){
			
			weight = 1;
			if(!SorB){
			
				event = (LHEFEvent*) branchEvent->At(0);
				weight = event->Weight;
			}
		
			tempSelection += weight;		
				
			// Choose the light jets to be the jets with highest pt
			Jet *tempJet1, *tempJet2, *tempJet3, *tempJet4, *tempJet5, *tempJet6;
				
			tempJet1 = Jet1;
			tempJet2 = Jet2;
			tempJet3 = Jet3;
			tempJet4 = Jet4;
			tempJet5 = Jet5;
			tempJet6 = Jet6;
			
			Jet1 = tempJet1;
			Jet2 = tempJet2;
			bJet1 = tempJet3;
			bJet2 = tempJet4;
			bJet3 = tempJet5;
			bJet4 = tempJet6;
			
			// Set up four-vectors for light jets
			vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
			vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
			
			vdiJet = vJet1 + vJet2;
			
			(SorB ? hMassdiJetS : hMassdiJetB)->Fill(vdiJet.M());
			
			// Check for dijet cut requirements
			if((vdiJet.M() > 500) && (fabs(Jet1->Eta - Jet2->Eta) > 2.5)){
				
				tempDijetCut += weight;
				
				// Set up four-vectors for bjets
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
				
				(SorB ? hInJetsS : hInJetsB)->Fill(inBJets);
				
				// Check if all bjets are between light jets
				if(inBJets==4){
					
					tempInjetCut += weight;
					
					v4bJet = vbJet1 + vbJet2 + vbJet3 + vbJet4;
				
					(SorB ? hMass4bJetsS : hMass4bJetsB)->Fill(v4bJet.M());
					
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
					(SorB ? hMass2bJetsNoCrossS : hMass2bJetsNoCrossB)->Fill(v2bJetNoCross.M());
					v2bJetNoCross = vbJet03 + vbJet04;
					(SorB ? hMass2bJetsNoCrossS : hMass2bJetsNoCrossB)->Fill(v2bJetNoCross.M());
					
					v2bJetCross = vbJet01 + vbJet03;
					(SorB ? hMass2bJetsCrossS : hMass2bJetsCrossB)->Fill(v2bJetCross.M());
					v2bJetCross = vbJet02 + vbJet04;
					(SorB ? hMass2bJetsCrossS : hMass2bJetsCrossB)->Fill(v2bJetCross.M());
				}
				
			}
			
		}
		
	} // End event loop
	
	// Update global variables
	(SorB ? selectionSignal : selectionBackground) += crossSection*3000*tempSelection/numberOfEntries;
	(SorB ? dijetCutSignal : dijetCutBackground) += crossSection*3000*tempDijetCut/numberOfEntries;
	(SorB ? injetCutSignal : injetCutBackground) += crossSection*3000*tempInjetCut/numberOfEntries;
	
	ifs.close();
}

/*
 * FUNCTION FOR SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 900, 600);
	
	gStyle->SetOptStat(kFALSE);

	// Save histograms
	histogram(hMassdiJetS, c1, "Mass diJet system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMassdiJetsS_2.jpg");
	histogram(hMass4bJetsS, c1, "Mass 4 bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass4bJetsS_2.jpg");
	histogram(hMass2bJetsNoCrossS, c1, "Mass 2 not crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsNoCrossS_2.jpg");
	histogram(hMass2bJetsCrossS, c1, "Mass 2 crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsCrossS_2.jpg");
	histogram(hInJetsS, c1, "Number of bjets between light jets", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramInBJetsS_2.jpg");
	
	histogram(hMassdiJetB, c1, "Mass diJet system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMassdiJetsB_2.jpg");
	histogram(hMass4bJetsB, c1, "Mass 4 bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass4bJetsB_2.jpg");
	histogram(hMass2bJetsNoCrossB, c1, "Mass 2 not crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsNoCrossB_2.jpg");
	histogram(hMass2bJetsCrossB, c1, "Mass 2 crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsCrossB_2.jpg");
	histogram(hInJetsB, c1, "Number of bjets between light jets", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramInBJetsB_2.jpg");
	
	histogram(hMassdiJetS, hMassdiJetB, c1, "Mass diJet system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMassdiJets_2.jpg");
	histogram(hMass4bJetsS, hMass4bJetsB, c1, "Mass 4 bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass4bJets_2.jpg");
	histogram(hMass2bJetsNoCrossS, hMass2bJetsNoCrossB, c1, "Mass 2 not crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsNoCross_2.jpg");
	histogram(hMass2bJetsCrossS, hMass2bJetsCrossB, c1, "Mass 2 crossed bjets system", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramMass2bJetsCross_2.jpg");
	histogram(hInJetsS, hInJetsB, c1, "Number of bjets between light jets", "Count", "/afs/cern.ch/user/a/ariostas/vbf-bbbb/Histograms/histogramInBJets_2.jpg");
	
	// Print event yields
	cout << "Total events: " << 1.8*3000*0.577*0.577 << endl;
	cout << "Events with selected jets: " << selectionSignal << endl;
	cout << "Events after dijet cut: " << dijetCutSignal << endl;
	cout << "Events after injet cut: " << injetCutSignal << endl;
	
	cout << "\nEvents with selected jets: " << selectionBackground << endl;
	cout << "Events after dijet cut: " << dijetCutBackground << endl;
	cout << "Events after injet cut: " << injetCutBackground << endl << endl;

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

