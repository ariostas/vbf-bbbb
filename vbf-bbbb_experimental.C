/*
 * 
 * This code is experimental, it runs over the  small ntuples generated with
 * the experimental selection. It's very messy and buggy, use with caution.
 * 
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
#include <TLegend.h>
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
TH1D *hMassdiJetS = new TH1D("MassdiJetS", "MassdiJetS", 200, 0, 4000);
TH1D *hInJetsS = new TH1D("InJetsS", "InJetsS", 7, -0.5, 6.5);
TH1D *hdRJJS = new TH1D("hdRJJS", "hdRJJS", 200, 0, 6);
TH1D *hdRBJS = new TH1D("hdRBJS", "hdRBJS", 200, 0, 9);
TH1D *hdES = new TH1D("hdES", "hdES", 200, -0.5, 10.5);

TH1D *hTestS = new TH1D("testS", "testS", 200, 0, 1000);

TH1D *hMassdiJetB = new TH1D("MassdiJetB", "MassdiJetB", 200, 0, 4000);
TH1D *hInJetsB = new TH1D("InJetsB", "InJetsB", 7, -0.5, 6.5);
TH1D *hdRJJB = new TH1D("hdRJJB", "hdRJJB", 200, 0, 6);
TH1D *hdRBJB = new TH1D("hdRBJB", "hdRBJB", 200, 0, 9);
TH1D *hdEB = new TH1D("hdEB", "hdEB", 200, -0.5, 10.5);

TH1D *hTestB = new TH1D("testB", "testB", 200, 0, 1000);

TH1D *hTest1S = new TH1D("test1S", "test1S", 7, -0.5, 6.5);						TH1D *hTest1B = new TH1D("test1B", "test1B", 7, -0.5, 6.5);
TH1D *hTest2S = new TH1D("test2S", "test2S", 11, -0.5, 10.5);					TH1D *hTest2B = new TH1D("test2B", "test2B", 11, -0.5, 10.5);
TH1D *hTest3S = new TH1D("test3S", "test3S", 11, -0.5, 10.5);						TH1D *hTest3B = new TH1D("test3B", "test3B", 11, -0.5, 10.5);
TH1D *hTest4S = new TH1D("test4S", "test4S", 250, 0, 600);						TH1D *hTest4B = new TH1D("test4B", "test4B", 250, 0, 600);
TH1D *hTest5S = new TH1D("test5S", "test5S", 6, -0.5, 5.5);						TH1D *hTest5B = new TH1D("test5B", "test5B", 6, -0.5, 5.5);

TH1D *hTest6S = new TH1D("test6S", "test6S", 100, 30, 150);						TH1D *hTest6B = new TH1D("test6B", "test6B", 100, 30, 150);
TH1D *hTest7S = new TH1D("test7S", "test7S", 100, 30, 150);						TH1D *hTest7B = new TH1D("test7B", "test7B", 100, 30, 150);
TH1D *hTest8S = new TH1D("test8S", "test8S", 100, 30, 150);						TH1D *hTest8B = new TH1D("test8B", "test8B", 100, 30, 150);
TH1D *hTest9S = new TH1D("test9S", "test9S", 100, 30, 150);						TH1D *hTest9B = new TH1D("test9B", "test9B", 100, 30, 150);
TH1D *hTest10S = new TH1D("test10S", "test10S", 100, 0, 350);					TH1D *hTest10B = new TH1D("test10B", "test10B", 100, 0, 350);
TH1D *hTest11S = new TH1D("test11S", "test11S", 100, 30, 150);					TH1D *hTest11B = new TH1D("test11B", "test11B", 100, 30, 150);

TH1D *hTest12S = new TH1D("test12S", "test12S", 200, 0, 250);					TH1D *hTest12B = new TH1D("test12B", "test12B", 200, 0, 250);
TH1D *hTest13S = new TH1D("test13S", "test13S", 100, 0, 300);					TH1D *hTest13B = new TH1D("test13B", "test13B", 100, 0, 300);
TH1D *hTest14S = new TH1D("test14S", "test14S", 100, -5, 5);					TH1D *hTest14B = new TH1D("test14B", "test14B", 100, -5, 5);

// Initialyze storage variables
Double_t totalSignal=0, selectionSignal = 0, kinCutSignal = 0, massCutSignal = 0;
Double_t totalBackground=0, selectionBackground = 0, kinCutBackground = 0, massCutBackground = 0;

Double_t ErrorSelectionSignal = 0, ErrorKinCutSignal = 0, ErrorMassCutSignal = 0;
Double_t ErrorSelectionBackground = 0, ErrorKinCutBackground = 0, ErrorMassCutBackground = 0;

/*
 * MAIN FUNCTION
 */

void vbf_bbbb_experimental(TString backgroundSample = "all"){
	
	BackgroundSample = backgroundSample;
	
	// Analyze signal
	analyze("HHToBBBB_14TeV", 0.669, Signal);
	//analyze("gfHHToBBBB_14TeV", 40*0.577*0.577, Signal);
	
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
		analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
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
		analyze("LLB-4p-400-900-v1510_14TEV", 0.22854*1000, Background);
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
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/" + inputfile + ".root";
	
	inputfile = "Reading " + inputfile + " events... ";
	
	inputfile.Resize(60);
	
	cout << inputfile;

	// Set up storage variables
	Double_t bjet1_Pt, bjet1_Eta, bjet1_Phi, bjet1_Mass;
	Double_t bjet2_Pt, bjet2_Eta, bjet2_Phi, bjet2_Mass;
	Double_t bjet3_Pt, bjet3_Eta, bjet3_Phi, bjet3_Mass;
	Double_t bjet4_Pt, bjet4_Eta, bjet4_Phi, bjet4_Mass;
	Double_t jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass;
	Double_t jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass;
	Double_t weight;
	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vDiJet, vDiHiggs, v2bJetNoCross, v2bJetCross, vbJet01, vbJet02, vbJet03, vbJet04;
	Int_t nLeptons=0, nLeptons01=0, nLeptons04=0, nJets=0, nBJets=0, nJetsPU=0, nBJetsPU=0, nDiJets=0;

	TFile* infile = new TFile(inputFile); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("weight",		&weight);
	intree->SetBranchAddress("bjet1_Pt",	&bjet1_Pt);
	intree->SetBranchAddress("bjet1_Eta",	&bjet1_Eta);
	intree->SetBranchAddress("bjet1_Phi",	&bjet1_Phi);
	intree->SetBranchAddress("bjet1_Mass",	&bjet1_Mass);
	intree->SetBranchAddress("bjet2_Pt",	&bjet2_Pt);
	intree->SetBranchAddress("bjet2_Eta",	&bjet2_Eta);
	intree->SetBranchAddress("bjet2_Phi",	&bjet2_Phi);
	intree->SetBranchAddress("bjet2_Mass",	&bjet2_Mass);
	intree->SetBranchAddress("bjet3_Pt",	&bjet3_Pt);
	intree->SetBranchAddress("bjet3_Eta",	&bjet3_Eta);
	intree->SetBranchAddress("bjet3_Phi",	&bjet3_Phi);
	intree->SetBranchAddress("bjet3_Mass",	&bjet3_Mass);
	intree->SetBranchAddress("bjet4_Pt",	&bjet4_Pt);
	intree->SetBranchAddress("bjet4_Eta",	&bjet4_Eta);
	intree->SetBranchAddress("bjet4_Phi",	&bjet4_Phi);
	intree->SetBranchAddress("bjet4_Mass",	&bjet4_Mass);
	intree->SetBranchAddress("jet1_Pt",		&jet1_Pt);
	intree->SetBranchAddress("jet1_Eta",	&jet1_Eta);
	intree->SetBranchAddress("jet1_Phi",	&jet1_Phi);
	intree->SetBranchAddress("jet1_Mass",	&jet1_Mass);
	intree->SetBranchAddress("jet2_Pt",		&jet2_Pt);
	intree->SetBranchAddress("jet2_Eta",	&jet2_Eta);
	intree->SetBranchAddress("jet2_Phi",	&jet2_Phi);
	intree->SetBranchAddress("jet2_Mass",	&jet2_Mass);
	intree->SetBranchAddress("nLeptons",    &nLeptons);
	intree->SetBranchAddress("nLeptons01",  &nLeptons01);
	intree->SetBranchAddress("nLeptons04",  &nLeptons04);
	intree->SetBranchAddress("nJets",     	&nJets);
	intree->SetBranchAddress("nBJets",     	&nBJets);
	intree->SetBranchAddress("nJetsPU",     &nJetsPU);
	intree->SetBranchAddress("nBJetsPU",    &nBJetsPU);
	intree->SetBranchAddress("nDiJets",     &nDiJets);
	
	// Set up temporal variables
	Double_t tempSelection=0, tempKinCut=0, tempMassCut=0;
	Double_t tempErrorSelection=0, tempErrorKinCut=0, tempErrorMassCut=0;

	for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // Event loop
		intree->GetEntry(iEntry);
		
		// Set up four-vectors for light jets
		vJet1.SetPtEtaPhiM(jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass);
		vJet2.SetPtEtaPhiM(jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass);
		
		vbJet1.SetPtEtaPhiM(bjet1_Pt, bjet1_Eta, bjet1_Phi, bjet1_Mass);
		vbJet2.SetPtEtaPhiM(bjet2_Pt, bjet2_Eta, bjet2_Phi, bjet2_Mass);
		vbJet3.SetPtEtaPhiM(bjet3_Pt, bjet3_Eta, bjet3_Phi, bjet3_Mass);
		vbJet4.SetPtEtaPhiM(bjet4_Pt, bjet4_Eta, bjet4_Phi, bjet4_Mass);
	
		vDiJet = vJet1 + vJet2;
		vDiHiggs = vbJet1 + vbJet2 + vbJet3 + vbJet4;
		
		Int_t inBJets = 0;
		if(jet1_Eta > jet2_Eta){
			
			if((bjet1_Eta < jet1_Eta) && (bjet1_Eta > jet2_Eta))inBJets++;
			if((bjet2_Eta < jet1_Eta) && (bjet2_Eta > jet2_Eta))inBJets++;
			if((bjet3_Eta < jet1_Eta) && (bjet3_Eta > jet2_Eta))inBJets++;
			if((bjet4_Eta < jet1_Eta) && (bjet4_Eta > jet2_Eta))inBJets++;
			
		}
		else{
			
			if((bjet1_Eta > jet1_Eta) && (bjet1_Eta < jet2_Eta))inBJets++;
			if((bjet2_Eta > jet1_Eta) && (bjet2_Eta < jet2_Eta))inBJets++;
			if((bjet3_Eta > jet1_Eta) && (bjet3_Eta < jet2_Eta))inBJets++;
			if((bjet4_Eta > jet1_Eta) && (bjet4_Eta < jet2_Eta))inBJets++;
			
		}
			
		tempSelection += weight;
		tempErrorSelection++;
		
		Double_t mindRb, mindRb1, mindRb2, mindRb3, mindRb4, mindRb5, mindRb6;
		
		mindRb1 = deltaR(bjet1_Eta, bjet2_Eta, bjet1_Phi, bjet2_Phi);
		mindRb2 = deltaR(bjet1_Eta, bjet3_Eta, bjet1_Phi, bjet3_Phi);
		mindRb3 = deltaR(bjet1_Eta, bjet4_Eta, bjet1_Phi, bjet4_Phi);
		mindRb4 = deltaR(bjet2_Eta, bjet3_Eta, bjet2_Phi, bjet3_Phi);
		mindRb5 = deltaR(bjet2_Eta, bjet4_Eta, bjet2_Phi, bjet4_Phi);
		mindRb6 = deltaR(bjet3_Eta, bjet4_Eta, bjet3_Phi, bjet4_Phi);
			
		mindRb = mindRb1;
		if(mindRb2 < mindRb) mindRb = mindRb2;
		if(mindRb3 < mindRb) mindRb = mindRb3;
		if(mindRb4 < mindRb) mindRb = mindRb4;
		if(mindRb5 < mindRb) mindRb = mindRb5;
		if(mindRb6 < mindRb) mindRb = mindRb6;
		
		Double_t mindRbj, mindRbj1, mindRbj2, mindRbj3, mindRbj4, mindRbj5, mindRbj6, mindRbj7, mindRbj8;
		
		mindRbj1 = deltaR(jet1_Eta, bjet1_Eta, jet1_Phi, bjet1_Phi);
		mindRbj2 = deltaR(jet1_Eta, bjet2_Eta, jet1_Phi, bjet2_Phi);
		mindRbj3 = deltaR(jet1_Eta, bjet3_Eta, jet1_Phi, bjet3_Phi);
		mindRbj4 = deltaR(jet1_Eta, bjet4_Eta, jet1_Phi, bjet4_Phi);
		mindRbj5 = deltaR(jet2_Eta, bjet1_Eta, jet2_Phi, bjet1_Phi);
		mindRbj6 = deltaR(jet2_Eta, bjet2_Eta, jet2_Phi, bjet2_Phi);
		mindRbj7 = deltaR(jet2_Eta, bjet3_Eta, jet2_Phi, bjet3_Phi);
		mindRbj8 = deltaR(jet2_Eta, bjet4_Eta, jet2_Phi, bjet4_Phi);
			
		mindRbj = mindRbj1;
		if(mindRbj2 < mindRbj) mindRbj = mindRbj2;
		if(mindRbj3 < mindRbj) mindRbj = mindRbj3;
		if(mindRbj4 < mindRbj) mindRbj = mindRbj4;
		if(mindRbj5 < mindRbj) mindRbj = mindRbj5;
		if(mindRbj6 < mindRbj) mindRbj = mindRbj6;
		if(mindRbj7 < mindRbj) mindRbj = mindRbj7;
		if(mindRbj8 < mindRbj) mindRbj = mindRbj8;
		
		Double_t dRJ = deltaR(jet1_Eta, jet2_Eta, jet1_Eta, jet2_Eta);
			
		tempKinCut += weight;
		tempErrorKinCut++;	
			
		// Set up temporal TLorentzVectors
		TLorentzVector test11, test12, test21, test22, test31, test32, higgs1, higgs2;
		test11 = vbJet1 + vbJet2;
		test12 = vbJet3 + vbJet4;
		test21 = vbJet1 + vbJet3;
		test22 = vbJet2 + vbJet4;
		test31 = vbJet1 + vbJet4;
		test32 = vbJet3 + vbJet2;
		
		Int_t nHiggs=0;
		
		if((test11.M() > 90) && (test11.M() < 135) && (test12.M() > 90) && (test12.M() < 135)) nHiggs++; 
		if((test21.M() > 90) && (test21.M() < 135) && (test22.M() > 90) && (test22.M() < 135)) nHiggs++;
		if((test31.M() > 90) && (test31.M() < 135) && (test32.M() > 90) && (test32.M() < 135)) nHiggs++;
		
		//if((test11.Pt() > 50) && (test12.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); (SorB ? hTest10S : hTest10B)->Fill(test11.Pt(), weight); (SorB ? hTest10S : hTest10B)->Fill(test12.Pt(), weight);}
		//if((test21.Pt() > 50) && (test22.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); (SorB ? hTest10S : hTest10B)->Fill(test21.Pt(), weight); (SorB ? hTest10S : hTest10B)->Fill(test22.Pt(), weight);}
		//if((test31.Pt() > 50) && (test32.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); (SorB ? hTest10S : hTest10B)->Fill(test31.Pt(), weight); (SorB ? hTest10S : hTest10B)->Fill(test32.Pt(), weight);}
		
		//if((test11.M() > 85) && (test11.M() < 145) && (test12.M() > 85) && (test12.M() < 145) && (test11.Pt() > 50) && (test12.Pt() > 50) && (mindRb1 > 0.5) && (mindRb1 < 2) && (mindRb6 > 0.5) && (mindRb6 < 2)){ nHiggs++; higgs1 = test11; higgs2 = test12;}
		//if((test21.M() > 85) && (test21.M() < 145) && (test22.M() > 85) && (test22.M() < 145) && (test21.Pt() > 50) && (test22.Pt() > 50) && (mindRb2 > 0.5) && (mindRb2 < 2) && (mindRb5 > 0.5) && (mindRb5 < 2)){ nHiggs++; higgs1 = test21; higgs2 = test22;}
		//if((test31.M() > 85) && (test31.M() < 145) && (test32.M() > 85) && (test32.M() < 145) && (test31.Pt() > 50) && (test32.Pt() > 50) && (mindRb3 > 0.5) && (mindRb3 < 2) && (mindRb4 > 0.5) && (mindRb4 < 2)){ nHiggs++; higgs1 = test31; higgs2 = test32;}
		
		
		//if((test11.M() > 90) && (test11.M() < 140) && (test12.M() > 90) && (test12.M() < 140) && (test11.Pt() > 50) && (test12.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); higgs1 = test11; higgs2 = test12;}
		//if((test21.M() > 90) && (test21.M() < 140) && (test22.M() > 90) && (test22.M() < 140) && (test21.Pt() > 50) && (test22.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); higgs1 = test21; higgs2 = test22;}
		//if((test31.M() > 90) && (test31.M() < 140) && (test32.M() > 90) && (test32.M() < 140) && (test31.Pt() > 50) && (test32.Pt() > 50)){ nHiggs++; (SorB ? hdRJJS : hdRJJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight); higgs1 = test31; higgs2 = test32;}
		
		Int_t ndRbb=0;
		
		//if((mindRb1 > 0.5) && (mindRb1 < 2) && (mindRb6 > 0.5) && (mindRb6 < 2)){ ndRbb++; (SorB ? hTest12S : hTest12B)->Fill(test11.M(), weight); (SorB ? hTest12S : hTest12B)->Fill(test12.M(), weight);}
		//if((mindRb2 > 0.5) && (mindRb2 < 2) && (mindRb5 > 0.5) && (mindRb5 < 2)){ ndRbb++; (SorB ? hTest12S : hTest12B)->Fill(test21.M(), weight); (SorB ? hTest12S : hTest12B)->Fill(test22.M(), weight);}
		//if((mindRb3 > 0.5) && (mindRb3 < 2) && (mindRb4 > 0.5) && (mindRb4 < 2)){ ndRbb++; (SorB ? hTest12S : hTest12B)->Fill(test31.M(), weight); (SorB ? hTest12S : hTest12B)->Fill(test32.M(), weight);}
		
		////////////////////////////////////////////////////////////////
		
		vJet1.SetPtEtaPhiM(jet1_Pt, jet1_Eta, jet1_Phi, jet1_Mass);
		vJet2.SetPtEtaPhiM(jet2_Pt, jet2_Eta, jet2_Phi, jet2_Mass);
		
		if(deltaR(jet1_Eta, test11.Eta(), jet1_Phi, test11.Phi()) < deltaR(jet2_Eta, test11.Eta(), jet2_Phi, test11.Phi())){
			
			test21 = test11 + vJet1;
			test22 = test12 + vJet2;
			
		}
		
		else{
			
			test21 = test11 + vJet2;
			test22 = test12 + vJet1;
		}
		
		if(nLeptons04 != 0) continue;
		if(nJetsPU > 7) continue;
		
		//if(bjet1_Pt < 30) continue;
		//if(bjet2_Pt < 30) continue;
		//if(bjet3_Pt < 30) continue;
		//if(bjet4_Pt < 30) continue;
		//if(jet1_Pt < 30) continue;
		//if(jet2_Pt < 30) continue;
		
		if(inBJets!=4) continue;
		
		if(mindRb < 0.75) continue;
		if(mindRbj < 1.5) continue;
		if(fabs(jet1_Eta - jet2_Eta) < 4) continue;
		//if(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()) < 2) continue;
		
		//if(vDiHiggs.M() < 275) continue;
		//if(nHiggs!=1) continue;
		//if(vDiJet.M() < 500) continue;
		
		//if(nCentralLightJets != 0) continue;
		
		//if(ndRbb!=1) continue;
		
		if(test11.M() <= 90 || test11.M() >= 150) continue;
		if(test12.M() <= 90 || test12.M() >= 150) continue;
		
		if(bjet1_Pt <= 80 || bjet3_Pt <= 80) continue;
		
		tempMassCut += weight;
		tempErrorMassCut++;
		
		
		(SorB ? hTest3S : hTest3B)->Fill(nJetsPU, weight);
		//(SorB ? hTest2S : hTest2B)->Fill(nCentralLightJets, weight);
		(SorB ? hTest5S : hTest5B)->Fill(nLeptons04, weight);
		
		(SorB ? hTest6S : hTest6B)->Fill(bjet1_Pt, weight);
		(SorB ? hTest7S : hTest7B)->Fill(bjet2_Pt, weight);
		(SorB ? hTest8S : hTest8B)->Fill(bjet3_Pt, weight);
		(SorB ? hTest9S : hTest9B)->Fill(bjet4_Pt, weight);
		(SorB ? hTest10S : hTest10B)->Fill(jet1_Pt, weight);
		(SorB ? hTest11S : hTest11B)->Fill(jet2_Pt, weight);
		(SorB ? hInJetsS : hInJetsB)->Fill(inBJets, weight);
		(SorB ? hTest12S : hTest12B)->Fill(test11.M(), weight);
		(SorB ? hTest12S : hTest12B)->Fill(test12.M(), weight);
		(SorB ? hTest13S : hTest13B)->Fill(vDiJet.Pt(), weight);
		(SorB ? hTest14S : hTest14B)->Fill(vDiHiggs.Eta(), weight);
		(SorB ? hMassdiJetS : hMassdiJetB)->Fill(vDiJet.M(), weight);
		(SorB ? hTest1S : hTest1B)->Fill(nHiggs, weight);
		(SorB ? hdRJJS : hdRJJB)->Fill(mindRbj, weight);
		(SorB ? hdRBJS : hdRBJB)->Fill(deltaR(test11.Eta(), test12.Eta(), test11.Phi(), test12.Phi()), weight);
		//(SorB ? hdRBJS : hdRBJB)->Fill(fabs(jet2_Eta), weight);
		(SorB ? hdES : hdEB)->Fill(fabs(jet1_Eta - jet2_Eta), weight);
		//(SorB ? hdES : hdEB)->Fill(fabs(test11.Phi() - test12.Phi()), weight);
		

		////////////////////////////////temp////////////////////////////
		
		/*if(nLeptons01 != 0) continue;
		if(nJets != 6) continue;
		if(bjet1_Pt < 40) continue;
		if(bjet2_Pt < 40) continue;
		if(bjet3_Pt < 30) continue;
		if(bjet4_Pt < 30) continue;
		if(jet1_Pt < 30) continue;
		if(jet2_Pt < 30) continue;
		if(inBJets!=4) continue;
		
		if(fabs(jet1_Eta - jet2_Eta) < 4.25) continue;
		if(mindRb < 0.75) continue;
		if(nHiggs!=1) continue;
		if(vDiHiggs.M() < 275) continue;
		if(vDiJet.M() < 400) continue;*/
		
		////////////////////////////////////////////////////////////////
		
	} // end event loop
	
	TString out = "";
	out += tempMassCut;
	out.Resize(8);
		
	cout << out << " events passed all cuts" << endl;
	
	// Update global variables
	(SorB ? totalSignal : totalBackground) += 3000*crossSection;
	(SorB ? selectionSignal : selectionBackground) += tempSelection;
	(SorB ? kinCutSignal : kinCutBackground) += tempKinCut;
	(SorB ? massCutSignal : massCutBackground) += tempMassCut;
	
	if(tempErrorSelection > 0) (SorB ? ErrorSelectionSignal : ErrorSelectionBackground) += sqrtf(tempErrorSelection)*tempSelection/tempErrorSelection;
	if(tempErrorKinCut > 0) (SorB ? ErrorKinCutSignal : ErrorKinCutBackground) += sqrtf(tempErrorKinCut)*tempKinCut/tempErrorKinCut;
	if(tempErrorMassCut > 0) (SorB ? ErrorMassCutSignal : ErrorMassCutBackground) += sqrtf(tempErrorMassCut)*tempMassCut/tempErrorMassCut;
	
}

/*
 * FUNCTION FOR SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);


	histogram(hdRJJS, hdRJJB, c1, "Minimum deltaR between bjets", "Count", "./Histograms/histogramdRbb_" + BackgroundSample + ".jpg");
	histogram(hdRBJS, hdRBJB, c1, "DeltaR between light jets", "Count", "./Histograms/histogramdRjj_" + BackgroundSample + ".jpg");
	histogram(hdES, hdEB, c1, "Delta eta between light jets", "Count", "./Histograms/histogramdEJ_" + BackgroundSample + ".jpg");
	histogram(hMassdiJetS, hMassdiJetB, c1, "Reconstructed DiJet mass", "Count", "./Histograms/histogramMassdiJets_" + BackgroundSample + ".jpg");
	histogram(hInJetsS, hInJetsB, c1, "Number of bjets between light jets", "Count", "./Histograms/histogramInBJets_" + BackgroundSample + ".jpg");
	
	histogram(hTest1S, hTest1B, c1, "Number of combinations with both M_bb in the range 90-135", "Count", "./Histograms/histogramTest1_" + BackgroundSample + ".jpg");
	histogram(hTest3S, hTest3B, c1, "Number of non-pileup jets with pt > 30", "Count", "./Histograms/histogramTest3_" + BackgroundSample + ".jpg");
	histogram(hTest5S, hTest5B, c1, "Number of leptons", "Count", "./Histograms/histogramTest5_" + BackgroundSample + ".jpg");
	
	histogram(hTest2S, hTest2B, c1, "Number of high pt jets", "Count", "./Histograms/histogramTest2_" + BackgroundSample + ".jpg");
	
	histogram(hTest6S, hTest6B, c1, "bJet1 pt", "Count", "./Histograms/histogramTest6_" + BackgroundSample + ".jpg");
	histogram(hTest7S, hTest7B, c1, "bJet2 pt", "Count", "./Histograms/histogramTest7_" + BackgroundSample + ".jpg");
	histogram(hTest8S, hTest8B, c1, "bJet3 pt", "Count", "./Histograms/histogramTest8_" + BackgroundSample + ".jpg");
	histogram(hTest9S, hTest9B, c1, "bJet4 pt", "Count", "./Histograms/histogramTest9_" + BackgroundSample + ".jpg");
	histogram(hTest10S, hTest10B, c1, "Jet1 pt", "Count", "./Histograms/histogramTest10_" + BackgroundSample + ".jpg");
	histogram(hTest11S, hTest11B, c1, "Jet2 pt", "Count", "./Histograms/histogramTest11_" + BackgroundSample + ".jpg");
	
	histogram(hTest12S, hTest12B, c1, "Reconstructed DiHiggs mass", "Count", "./Histograms/histogramTest12_" + BackgroundSample + ".jpg");
	histogram(hTest13S, hTest13B, c1, "HH pt", "Count", "./Histograms/histogramTest13_" + BackgroundSample + ".jpg");
	histogram(hTest14S, hTest14B, c1, "HH Eta", "Count", "./Histograms/histogramTest14_" + BackgroundSample + ".jpg");
	
	// Print event yields
	cout << "\nSignal" << endl << endl;
	cout << "Total events: " << totalSignal << endl;
	cout << "Events after selection: " << selectionSignal << " +- " << ErrorSelectionSignal << endl;
	cout << "Events after kinematic cuts: " << kinCutSignal << " +- " << ErrorKinCutSignal << endl;
	cout << "Events after mass cuts: " << massCutSignal << " +- " << ErrorMassCutSignal << endl;
	
	cout << BackgroundSample << " background" << endl << endl;
	cout << "Total events: " << totalBackground << endl;
	cout << "Events after selection: " << selectionBackground << " +- " << ErrorSelectionBackground << endl;
	cout << "Events after kinematic cuts: " << kinCutBackground << " +- " << ErrorKinCutBackground << endl;
	cout << "Events after mass cuts: " << massCutBackground << " +- " << ErrorMassCutBackground << endl;

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
	cout << "Events after event selection: " << selectionSignal << " +- " << ErrorSelectionSignal << endl;
	cout << "Events after kinematic cuts: " << kinCutSignal << " +- " << ErrorKinCutSignal << endl;
	cout << "Events after mass cuts: " << massCutSignal << " +- " << ErrorMassCutSignal << endl;

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
	
	TLegend *leg = new TLegend(0.6,0.65,0.88,0.85);
	leg->SetTextFont(72);
	leg->SetTextSize(0.04);
	leg->AddEntry(histoS,"Signal","l");
	leg->AddEntry(histoB, (BackgroundSample == "all" ?  "All backgrounds" : BackgroundSample + " background"),"l");
	leg->Draw();

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
