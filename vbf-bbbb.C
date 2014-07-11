/*
 * vbf-bbbb.C
 * 
 * This macro analyses signal and background samples for vbf-bbbb.
 * It analyses the samples by selecting 4 bjets and 2 light jets.
 * To run use "root vbf-bbbb_0.C+\(\"name of sample\"\)"
 * 
 * The selection of events consist in requiring all jets to pass 
 * the jet id veto, have pt > 40 and |eta| < 5
 * 
 * The dijet cut consist in requiring the light jets to yield a reconstructed
 * mass > 800 GeV and to have delta eta > 3.5
 * 
 * The injet cut requires all bjets to be between the light jets
 * 
 * The Higgs cut requires that at least one combination of pairs of bjets
 * yield a reconstructed mass in the range 110-140.
 * 
 * The dRb cut requires the minimim delta R between bjets to be greater than 1.
 * 
 * The dRj cut requires the delta R between light jets to be greater than 6.
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

TString BackgroundSample;

// Declare functions
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void analyze(TString, Double_t, bool);
void saveResults();
void saveResultsS();
int puJetID(Float_t, Float_t, Float_t);
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

// Initialize histograms
TH1D *hMassdiJetS = new TH1D("MassdiJetS", "MassdiJetS", 200, 0, 2600);
TH1D *hInJetsS = new TH1D("InJetsS", "InJetsS", 6, -0.5, 5.5);
TH1D *hdRJJS = new TH1D("hdRJJS", "hdRJJS", 100, -0.5, 10.5);
TH1D *hdRBJS = new TH1D("hdRBJS", "hdRBJS", 100, -0.5, 10.5);
TH1D *hdES = new TH1D("hdES", "hdES", 100, -0.5, 10.5);

TH1D *hTestS = new TH1D("testS", "testS", 200, 0, 1000);

TH1D *hMassdiJetB = new TH1D("MassdiJetB", "MassdiJetB", 200, 0, 2600);
TH1D *hInJetsB = new TH1D("InJetsB", "InJetsB", 6, -0.5, 5.5);
TH1D *hdRJJB = new TH1D("hdRJJB", "hdRJJB", 100, -0.5, 10.5);
TH1D *hdRBJB = new TH1D("hdRBJB", "hdRBJB", 100, -0.5, 10.5);
TH1D *hdEB = new TH1D("hdEB", "hdEB", 100, -0.5, 10.5);

TH1D *hTestB = new TH1D("testB", "testB", 200, 0, 1000);

// Initialyze storage variables
Double_t totalSignal=0, selectionSignal = 0, dijetCutSignal = 0, injetCutSignal = 0, higgsCutSignal=0, dRbCutSignal=0, dRjCutSignal=0;
Double_t totalBackground=0, selectionBackground = 0, dijetCutBackground = 0, injetCutBackground = 0, higgsCutBackground=0, dRbCutBackground=0, dRjCutBackground=0;

Double_t ErrorselectionSignal = 0, ErrordijetCutSignal = 0, ErrorinjetCutSignal = 0, ErrorhiggsCutSignal=0, ErrordRbCutSignal=0, ErrordRjCutSignal=0;
Double_t ErrorselectionBackground = 0, ErrordijetCutBackground = 0, ErrorinjetCutBackground = 0, ErrorhiggsCutBackground=0, ErrordRbCutBackground=0, ErrordRjCutBackground=0;

/*
 * MAIN FUNCTION
 */

void vbf_bbbb(TString backgroundSample){
	
	BackgroundSample = backgroundSample;
	
	// Analyze signal
	analyze("HHToBBBB_14TeV", 0.599, Signal);
	
	if(BackgroundSample == "B"){
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background);
	}
	
	else if(BackgroundSample == "BB"){
		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background);
	}
	
	else if(BackgroundSample == "BBB"){
		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background);
	}
	
	else if(BackgroundSample == "Bj"){
		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background);
	}
	
	else if(BackgroundSample == "Bjj"){
		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background);
	}
	
	else if(BackgroundSample == "H"){
		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background);
	}
	
	else if(BackgroundSample == "LL"){
		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background);
	}
	
	else if(BackgroundSample == "LLB"){
		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background);
	}
	
	else if(BackgroundSample == "tB"){
		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background);
	}
	
	else if(BackgroundSample == "tj"){
		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background);
	}
	
	else if(BackgroundSample == "tt"){
		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background);
	}
	
	else if(BackgroundSample == "ttB"){
		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background);
	}
	
	else if(BackgroundSample == "Signal"){
		cout << "Processing only signal events" << endl;
	}
	
	else if(BackgroundSample == "all"){
		
		// Analyze B background
		analyze("B-4p-0-1-v1510_14TEV", 200944.68129*1000, Background);

		// Analyse BB background
		analyze("BB-4p-0-300-v1510_14TEV", 249.97710*1000, Background);
		analyze("BB-4p-300-700-v1510_14TEV", 35.23062*1000, Background);
		analyze("BB-4p-700-1300-v1510_14TEV", 4.13743*1000, Background);
		analyze("BB-4p-1300-2100-v1510_14TEV", 0.41702*1000, Background);
		analyze("BB-4p-2100-100000-v1510_14TEV", 0.04770*1000, Background);

		// Analyze BBB background
		analyze("BBB-4p-0-600-v1510_14TEV", 2.57304*1000, Background);
		analyze("BBB-4p-600-1300-v1510_14TEV", 0.14935*1000, Background);
		analyze("BBB-4p-1300-100000-v1510_14TEV", 0.01274*1000, Background);

		// Analyze Bj background
		analyze("Bj-4p-0-300-v1510_14TEV", 34409.92339*1000, Background);
		analyze("Bj-4p-300-600-v1510_14TEV", 2642.85309*1000, Background);
		analyze("Bj-4p-600-1100-v1510_14TEV", 294.12311*1000, Background);
		analyze("Bj-4p-1100-1800-v1510_14TEV", 25.95000*1000, Background);
		analyze("Bj-4p-1800-2700-v1510_14TEV", 2.42111*1000, Background);
		analyze("Bj-4p-2700-3700-v1510_14TEV", 0.22690*1000, Background);
		analyze("Bj-4p-3700-100000-v1510_14TEV", 0.02767*1000, Background);

		// Analyze Bjj-vbf background
		analyze("Bjj-vbf-4p-0-700-v1510_14TEV", 86.45604*1000, Background);
		analyze("Bjj-vbf-4p-700-1400-v1510_14TEV", 4.34869*1000, Background);
		analyze("Bjj-vbf-4p-1400-2300-v1510_14TEV", 0.32465*1000, Background);
		analyze("Bjj-vbf-4p-2300-3400-v1510_14TEV", 0.03032*1000, Background);
		//analyze("Bjj-vbf-4p-3400-100000-v1510_14TEV", 0.00313*1000, Background);

		// Analyze H background
		analyze("H-4p-0-300-v1510_14TEV", 21.55990*1000, Background);
		analyze("H-4p-300-800-v1510_14TEV", 1.11282*1000, Background);
		analyze("H-4p-800-1500-v1510_14TEV", 0.09188*1000, Background);
		analyze("H-4p-1500-100000-v1510_14TEV", 0.01009*1000, Background);

		// Analyze LL background
		analyze("LL-4p-0-100-v1510_14TEV", 1341.36923*1000, Background);
		analyze("LL-4p-100-200-v1510_14TEV", 156.29534*1000, Background);
		analyze("LL-4p-200-500-v1510_14TEV", 42.40132*1000, Background);
		analyze("LL-4p-500-900-v1510_14TEV", 2.84373*1000, Background);
		analyze("LL-4p-900-1400-v1510_14TEV", 0.20914*1000, Background);
		analyze("LL-4p-1400-100000-v1510_14TEV", 0.02891*1000, Background);

		// Analyze LLB background
		analyze("LLB-4p-0-400-v1510_14TEV", 2.97380*1000, Background);
		//analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
		analyze("LLB-4p-900-100000-v1510_14TEV", 0.02080*1000, Background);

		// Analyze tB background
		analyze("tB-4p-0-500-v1510_14TEV", 63.88923*1000, Background);
		analyze("tB-4p-500-900-v1510_14TEV", 7.12172*1000, Background);
		analyze("tB-4p-900-1500-v1510_14TEV", 0.98030*1000, Background);
		analyze("tB-4p-1500-2200-v1510_14TEV", 0.08391*1000, Background);
		analyze("tB-4p-2200-100000-v1510_14TEV", 0.00953*1000, Background);

		// Analyze tj background
		analyze("tj-4p-0-500-v1510_14TEV", 109.73602*1000, Background);
		analyze("tj-4p-500-1000-v1510_14TEV", 5.99325*1000, Background);
		analyze("tj-4p-1000-1600-v1510_14TEV", 0.37680*1000, Background);
		analyze("tj-4p-1600-2400-v1510_14TEV", 0.03462*1000, Background);
		analyze("tj-4p-2400-100000-v1510_14TEV", 0.00312*1000, Background);

		// Analyse tt background
		analyze("tt-4p-0-600-v1510_14TEV", 530.89358*1000, Background);
		analyze("tt-4p-600-1100-v1510_14TEV", 42.55351*1000, Background);
		analyze("tt-4p-1100-1700-v1510_14TEV", 4.48209*1000, Background);
		analyze("tt-4p-1700-2500-v1510_14TEV", 0.52795*1000, Background);
		analyze("tt-4p-2500-100000-v1510_14TEV", 0.05449*1000, Background);

		// Analyze ttB background
		analyze("ttB-4p-0-900-v1510_14TEV", 2.6673*1000, Background);
		analyze("ttB-4p-900-1600-v1510_14TEV", 0.250469*1000, Background);
		analyze("ttB-4p-1600-2500-v1510_14TEV", 0.0237441*1000, Background);
		analyze("ttB-4p-2500-100000-v1510_14TEV", 0.00208816*1000, Background);
		
	}
	
	else {
		cout << "Background sample not found" << endl;
		assert(false);
	}
	
	
	// Save results
	if(BackgroundSample == "Signal") saveResultsS();
	else saveResults();
	
}

/*
 * ANALYSIS
 */

void analyze(TString inputfile, Double_t crossSection, bool SorB)
{	
	// The SorB variable is true if it's signal and false if it's a background
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/vbf-bbbb/" + inputfile + ".root";
	
	cout << "Reading from " << inputfile << endl;

	// Set up storage variables
	Jet *bJet1=0, *bJet2=0, *bJet3=0, *bJet4=0, *Jet1=0, *Jet2=0;
	LHEFEvent *event;
	Double_t weight;
	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vdiJet, v4bJet, v2bJetNoCross, v2bJetCross, vbJet01, vbJet02, vbJet03, vbJet04;

	TFile* infile = new TFile(inputFile); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("weight",	&weight);
	intree->SetBranchAddress("bJet1",		&bJet1);
	intree->SetBranchAddress("bJet2",		&bJet2);
	intree->SetBranchAddress("bJet3",		&bJet3);
	intree->SetBranchAddress("bJet4",		&bJet4);
	intree->SetBranchAddress("Jet1",     &Jet1);
	intree->SetBranchAddress("Jet2",     &Jet2);

	// Set up temporal variables
	Double_t tempSelection=0, tempDijetCut=0, tempInjetCut=0, tempHiggsCut=0, tempDRbCut=0, tempDRjCut=0;
	Double_t tempErrorSelection=0, tempErrorDijetCut=0, tempErrorInjetCut=0, tempErrorHiggsCut=0, tempErrorDRbCut=0, tempErrorDRjCut=0;

	for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // Event loop
		intree->GetEntry(iEntry);
			
		tempSelection += weight;
		tempErrorSelection++;
			
		// Set up four-vectors for light jets
		vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
		vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
		
		vdiJet = vJet1 + vJet2;
		
		(SorB ? hMassdiJetS : hMassdiJetB)->Fill(vdiJet.M(), weight);
		(SorB ? hdES : hdEB)->Fill(fabs(Jet1->Eta - Jet2->Eta), weight);
		
		// Check if it satisfies the dijet cut requirements
		if(!((vdiJet.M() > 800) && (fabs(Jet1->Eta - Jet2->Eta) > 3.5))) continue;
			
		tempDijetCut += weight;
		tempErrorDijetCut++;
		
		vbJet1.SetPtEtaPhiM(bJet1->PT, bJet1->Eta, bJet1->Phi, bJet1->Mass);
		vbJet2.SetPtEtaPhiM(bJet2->PT, bJet2->Eta, bJet2->Phi, bJet2->Mass);
		vbJet3.SetPtEtaPhiM(bJet3->PT, bJet3->Eta, bJet3->Phi, bJet3->Mass);
		vbJet4.SetPtEtaPhiM(bJet4->PT, bJet4->Eta, bJet4->Phi, bJet4->Mass);
	
		v4bJet = vbJet1 + vbJet2 + vbJet3 + vbJet4;
		
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
		
		(SorB ? hInJetsS : hInJetsB)->Fill(inBJets, weight);
		
		// Check if all bjets are between light jets
		if(!(inBJets==4)) continue;
		
		tempInjetCut += weight;	
		tempErrorInjetCut++;	
			
		TLorentzVector test11, test12, test21, test22, test31, test32;
		test11 = vbJet1 + vbJet2;
		test12 = vbJet3 + vbJet4;
		test21 = vbJet1 + vbJet3;
		test22 = vbJet2 + vbJet4;
		test31 = vbJet1 + vbJet4;
		test32 = vbJet3 + vbJet2;
		
		Int_t nHiggs=0;
		
		if((test11.M() > 110) && (test11.M() < 140) && (test12.M() > 110) && (test12.M() < 140)) nHiggs++; 
		if((test21.M() > 110) && (test21.M() < 140) && (test22.M() > 110) && (test22.M() < 140)) nHiggs++;
		if((test31.M() > 110) && (test31.M() < 140) && (test32.M() > 110) && (test32.M() < 140)) nHiggs++;
			
		//Check if at least one combination of pairs of bjets yield a reconstructed mass in the range 110-140
		if(!(nHiggs>0)) continue;
			
		tempHiggsCut += weight;
		tempErrorHiggsCut++;
			 
		Double_t mindRb, mindRb1, mindRb2, mindRb3, mindRb4, mindRb5, mindRb6;
		
		mindRb1 = deltaR(bJet1->Eta, bJet2->Eta, bJet1->Phi, bJet2->Phi);
		mindRb2 = deltaR(bJet1->Eta, bJet3->Eta, bJet1->Phi, bJet3->Phi);
		mindRb3 = deltaR(bJet1->Eta, bJet4->Eta, bJet1->Phi, bJet4->Phi);
		mindRb4 = deltaR(bJet2->Eta, bJet3->Eta, bJet2->Phi, bJet3->Phi);
		mindRb5 = deltaR(bJet2->Eta, bJet4->Eta, bJet2->Phi, bJet4->Phi);
		mindRb6 = deltaR(bJet3->Eta, bJet4->Eta, bJet3->Phi, bJet4->Phi);
			
		mindRb = mindRb1;
		if(mindRb2 < mindRb) mindRb = mindRb2;
		if(mindRb3 < mindRb) mindRb = mindRb3;
		if(mindRb4 < mindRb) mindRb = mindRb4;
		if(mindRb5 < mindRb) mindRb = mindRb5;
		if(mindRb6 < mindRb) mindRb = mindRb6;

		(SorB ? hdRJJS : hdRJJB)->Fill(mindRb, weight);
	
		//Check if it satisfies the min delta R_bb requirement
		if(!(mindRb > 1)) continue;
		
		tempDRbCut += weight;
		tempErrorDRbCut++;
					
		Double_t dRJ = deltaR(Jet1->Eta, Jet2->Eta, Jet1->Eta, Jet2->Eta);
					
		(SorB ? hdRBJS : hdRBJB)->Fill(dRJ, weight);
					
		if(!(dRJ > 6)) continue;
			
		tempDRjCut +=weight;
		tempErrorDRjCut++;		
		
	} // end event loop
	
	// Update global variables
	(SorB ? totalSignal : totalBackground) += 3000*crossSection;
	(SorB ? selectionSignal : selectionBackground) += tempSelection;
	(SorB ? dijetCutSignal : dijetCutBackground) += tempDijetCut;
	(SorB ? injetCutSignal : injetCutBackground) += tempInjetCut;
	(SorB ? higgsCutSignal : higgsCutBackground) += tempHiggsCut;
	(SorB ? dRbCutSignal : dRbCutBackground) += tempDRbCut;
	(SorB ? dRjCutSignal : dRjCutBackground) += tempDRjCut;
	
	(SorB ? ErrorselectionSignal : ErrorselectionBackground) += sqrtf(tempErrorSelection)*weight;
	(SorB ? ErrordijetCutSignal : ErrordijetCutBackground) += sqrtf(tempErrorDijetCut)*weight;
	(SorB ? ErrorinjetCutSignal : ErrorinjetCutBackground) += sqrtf(tempErrorInjetCut)*weight;
	(SorB ? ErrorhiggsCutSignal : ErrorhiggsCutBackground) += sqrtf(tempErrorHiggsCut)*weight;
	(SorB ? ErrordRbCutSignal : ErrordRbCutBackground) += sqrtf(tempErrorDRbCut)*weight;
	(SorB ? ErrordRjCutSignal : ErrordRjCutBackground) += sqrtf(tempErrorDRjCut)*weight;
	
}

/*
 * FUNCTION FOR SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);


	histogram(hdRJJS, hdRJJB, c1, "minimum deltaR between bjets", "Count", "./Histograms/histogramdRbb_" + BackgroundSample + ".jpg");
	histogram(hdRBJS, hdRBJB, c1, "deltaR betwenn light jets", "Count", "./Histograms/histogramdRjj_" + BackgroundSample + ".jpg");
	histogram(hdES, hdEB, c1, "|d_eta| between light jets", "Count", "./Histograms/histogramdEJ_" + BackgroundSample + ".jpg");
	histogram(hMassdiJetS, hMassdiJetB, c1, "Mass diJet system", "Count", "./Histograms/histogramMassdiJets_" + BackgroundSample + ".jpg");
	histogram(hInJetsS, hInJetsB, c1, "Number of bjets between light jets", "Count", "./Histograms/histogramInBJets_" + BackgroundSample + ".jpg");
	
	// Print event yields
	cout << "\nSignal" << endl << endl;
	cout << "Total events: " << totalSignal << endl;
	cout << "Events with selected jets: " << selectionSignal << " +- " << ErrorselectionSignal << endl;
	cout << "Events after dijet cut: " << dijetCutSignal << " +- " << ErrordijetCutSignal << endl;
	cout << "Events after injet cut: " << injetCutSignal << " +- " << ErrorinjetCutSignal << endl;
	cout << "Events after higgs cut: " << higgsCutSignal << " +- " << ErrorhiggsCutSignal << endl;
	cout << "Events after dRb cut: " << dRbCutSignal << " +- " << ErrordRbCutSignal << endl;
	cout << "Events after dRj cut: " << dRjCutSignal << " +- " << ErrordRjCutSignal << endl << endl << endl;
	
	cout << BackgroundSample << " background" << endl << endl;
	cout << "Total events: " << totalBackground << endl;
	cout << "Events with selected jets: " << selectionBackground << " +- " << ErrorselectionBackground << endl;
	cout << "Events after dijet cut: " << dijetCutBackground << " +- " << ErrordijetCutBackground << endl;
	cout << "Events after injet cut: " << injetCutBackground << " +- " << ErrorinjetCutBackground << endl;
	cout << "Events after higgs cut: " << higgsCutBackground << " +- " << ErrorhiggsCutBackground << endl;
	cout << "Events after dRb cut: " << dRbCutBackground << " +- " << ErrordRbCutBackground << endl;
	cout << "Events after dRj cut: " << dRjCutBackground << " +- " << ErrordRjCutBackground << endl << endl;

}

void saveResultsS(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);

	histogram(hdRJJS, c1, "minimum deltaR between bjets", "Count", "./Histograms/histogramdRbb_" + BackgroundSample + ".jpg");
	histogram(hdRBJS, c1, "deltaR betwenn light jets", "Count", "./Histograms/histogramdRjj_" + BackgroundSample + ".jpg");
	histogram(hdES, c1, "|d_eta| between light jets", "Count", "./Histograms/histogramdEJ_" + BackgroundSample + ".jpg");
	histogram(hMassdiJetS, c1, "Mass diJet system", "Count", "./Histograms/histogramMassdiJets_" + BackgroundSample + ".jpg");
	histogram(hInJetsS, c1, "Number of bjets between light jets", "Count", "./Histograms/histogramInBJets_" + BackgroundSample + ".jpg");
	
	// Print event yields
	cout << "Signal" << endl << endl;
	cout << "Total events: " << 1.8*3000*0.577*0.577 << endl;
	cout << "Events with selected jets: " << selectionSignal << " +- " << ErrorselectionSignal << endl;
	cout << "Events after dijet cut: " << dijetCutSignal << " +- " << ErrordijetCutSignal << endl;
	cout << "Events after injet cut: " << injetCutSignal << " +- " << ErrorinjetCutSignal << endl;
	cout << "Events after higgs cut: " << higgsCutSignal << " +- " << ErrorhiggsCutSignal << endl;
	cout << "Events after dRb cut: " << dRbCutSignal << " +- " << ErrordRbCutSignal << endl;
	cout << "Events after dRj cut: " << dRjCutSignal << " +- " << ErrordRjCutSignal << endl << endl << endl;

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
	if((histoS->GetMaximum()) > (histoB->GetMaximum())) max=1.1*(histoS->GetMaximum());
	else max=1.1*(histoB->GetMaximum());
	
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

Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
