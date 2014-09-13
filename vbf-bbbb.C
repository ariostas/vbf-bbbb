/*
 * vbf-bbbb.C
 * 
 * This macro analyses signal and background samples for vbf-bbbb.
 * To run use "root vbf-bbbb.C+" or "root vbf-bbbb.C+\(\"name of sample\"\)"
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
TH1D *hNLeptonsS = new TH1D("hNLeptonsS", "hNLeptonsS", 6, -0.5, 5.5);			TH1D *hNLeptonsB = new TH1D("hNLeptonsB", "hNLeptonsB", 6, -0.5, 5.5);
TH1D *hNJetsPUS = new TH1D("hNJetsPUS", "hNJetsPUS", 6, 5.5, 11.5);			TH1D *hNJetsPUB = new TH1D("hNJetsPUB", "hNJetsPUB", 6, 5.5, 11.5);
TH1D *hNBJetsS = new TH1D("hNBJetsS", "hNBJetsS", 7, -0.5, 6.5);				TH1D *hNBJetsB = new TH1D("hNBJetsB", "hNBJetsB", 7, -0.5, 6.5);

TH1D *hDRbbS = new TH1D("hDRbbS", "hDRbbS", 50, 0, 3);							TH1D *hDRbbB = new TH1D("hDRbbB", "hDRbbB", 50, 0, 3);
TH1D *hDRbjS = new TH1D("hDRbjS", "hDRbjS", 50, 0, 3);							TH1D *hDRbjB = new TH1D("hDRbjB", "hDRbjB", 50, 0, 3);
TH1D *hDRjjS = new TH1D("hDRjjS", "hDRjjS", 50, 0, 9);							TH1D *hDRjjB = new TH1D("hDRjjB", "hDRjjB", 50, 0, 9);
TH1D *hDEjjS = new TH1D("hDEjjS", "hDEjjS", 50, 0, 10);							TH1D *hDEjjB = new TH1D("hDEjjB", "hDEjjB", 50, 0, 10);

TH1D *hMjjS = new TH1D("hMjjS", "hMjjS", 50, 0, 4000);							TH1D *hMjjB = new TH1D("hMjjB", "hMjjB", 50, 0, 4000);
TH1D *hMHHS = new TH1D("hMHHS", "hMHHS", 50, 100, 1000);						TH1D *hMHHB = new TH1D("hMHHB", "hMHHB", 50, 100, 1000);
TH1D *hMbbS = new TH1D("hMbbS", "hMbbS", 7, -0.5, 6.5);							TH1D *hMbbB = new TH1D("hMbbB", "hMbbB", 7, -0.5, 6.5);

// Initialyze storage variables
Double_t totalSignal=0, selectionSignal = 0, kinCutSignal = 0, massCutSignal = 0;
Double_t totalBackground=0, selectionBackground = 0, kinCutBackground = 0, massCutBackground = 0;

Double_t ErrorSelectionSignal = 0, ErrorKinCutSignal = 0, ErrorMassCutSignal = 0;
Double_t ErrorSelectionBackground = 0, ErrorKinCutBackground = 0, ErrorMassCutBackground = 0;

/*
 * MAIN FUNCTION
 */

void vbf_bbbb(TString backgroundSample = "all"){
	
	BackgroundSample = backgroundSample;
	
	// Analyze signal
	analyze("HHToBBBB_14TeV", 2.01*0.577*0.577, Signal);
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
	
	const TString inputFile = "/afs/cern.ch/work/a/ariostas/public/vbf-bbbb/" + inputfile + ".root";
	
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
	Int_t nLeptons=0, nLeptons01=0, nLeptons04=0, nJets=0, nJetsHighPt=0;

	TFile* infile = new TFile(inputFile); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("weight",			&weight);
	intree->SetBranchAddress("bjet1_Pt",		&bjet1_Pt);
	intree->SetBranchAddress("bjet1_Eta",		&bjet1_Eta);
	intree->SetBranchAddress("bjet1_Phi",		&bjet1_Phi);
	intree->SetBranchAddress("bjet1_Mass",		&bjet1_Mass);
	intree->SetBranchAddress("bjet2_Pt",		&bjet2_Pt);
	intree->SetBranchAddress("bjet2_Eta",		&bjet2_Eta);
	intree->SetBranchAddress("bjet2_Phi",		&bjet2_Phi);
	intree->SetBranchAddress("bjet2_Mass",		&bjet2_Mass);
	intree->SetBranchAddress("bjet3_Pt",		&bjet3_Pt);
	intree->SetBranchAddress("bjet3_Eta",		&bjet3_Eta);
	intree->SetBranchAddress("bjet3_Phi",		&bjet3_Phi);
	intree->SetBranchAddress("bjet3_Mass",		&bjet3_Mass);
	intree->SetBranchAddress("bjet4_Pt",		&bjet4_Pt);
	intree->SetBranchAddress("bjet4_Eta",		&bjet4_Eta);
	intree->SetBranchAddress("bjet4_Phi",		&bjet4_Phi);
	intree->SetBranchAddress("bjet4_Mass",		&bjet4_Mass);
	intree->SetBranchAddress("jet1_Pt",			&jet1_Pt);
	intree->SetBranchAddress("jet1_Eta",		&jet1_Eta);
	intree->SetBranchAddress("jet1_Phi",		&jet1_Phi);
	intree->SetBranchAddress("jet1_Mass",		&jet1_Mass);
	intree->SetBranchAddress("jet2_Pt",			&jet2_Pt);
	intree->SetBranchAddress("jet2_Eta",		&jet2_Eta);
	intree->SetBranchAddress("jet2_Phi",		&jet2_Phi);
	intree->SetBranchAddress("jet2_Mass",		&jet2_Mass);
	intree->SetBranchAddress("nLeptons",     	&nLeptons);
	intree->SetBranchAddress("nLeptons01",  	&nLeptons01);
	intree->SetBranchAddress("nLeptons04",    	&nLeptons04);
	intree->SetBranchAddress("nJets",     		&nJets);
	intree->SetBranchAddress("nJetsHighPt",   	&nJetsHighPt);

	// Set up temporary variables
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
		
		Int_t nBjetsBetweenLjets = 0;
		if(jet1_Eta > jet2_Eta){
			
			if((bjet1_Eta < jet1_Eta) && (bjet1_Eta > jet2_Eta))nBjetsBetweenLjets++;
			if((bjet2_Eta < jet1_Eta) && (bjet2_Eta > jet2_Eta))nBjetsBetweenLjets++;
			if((bjet3_Eta < jet1_Eta) && (bjet3_Eta > jet2_Eta))nBjetsBetweenLjets++;
			if((bjet4_Eta < jet1_Eta) && (bjet4_Eta > jet2_Eta))nBjetsBetweenLjets++;
			
		}
		else{
			
			if((bjet1_Eta > jet1_Eta) && (bjet1_Eta < jet2_Eta))nBjetsBetweenLjets++;
			if((bjet2_Eta > jet1_Eta) && (bjet2_Eta < jet2_Eta))nBjetsBetweenLjets++;
			if((bjet3_Eta > jet1_Eta) && (bjet3_Eta < jet2_Eta))nBjetsBetweenLjets++;
			if((bjet4_Eta > jet1_Eta) && (bjet4_Eta < jet2_Eta))nBjetsBetweenLjets++;
			
		}
		
		(SorB ? hNBJetsS : hNBJetsB)->Fill(nBjetsBetweenLjets, weight);
		(SorB ? hNJetsPUS : hNJetsPUB)->Fill(nJets, weight);
		(SorB ? hNLeptonsS : hNLeptonsB)->Fill(nLeptons04, weight);
		
		//Check selection requirements
		if(nLeptons04 != 0) continue;
		if(nJets != 6) continue;
		if(nBjetsBetweenLjets!=4) continue;
			
		tempSelection += weight;
		tempErrorSelection++;
		
		Double_t minDRbb, minDRbb1, minDRbb2, minDRbb3, minDRbb4, minDRbb5, minDRbb6;
		
		minDRbb1 = deltaR(bjet1_Eta, bjet2_Eta, bjet1_Phi, bjet2_Phi);
		minDRbb2 = deltaR(bjet1_Eta, bjet3_Eta, bjet1_Phi, bjet3_Phi);
		minDRbb3 = deltaR(bjet1_Eta, bjet4_Eta, bjet1_Phi, bjet4_Phi);
		minDRbb4 = deltaR(bjet2_Eta, bjet3_Eta, bjet2_Phi, bjet3_Phi);
		minDRbb5 = deltaR(bjet2_Eta, bjet4_Eta, bjet2_Phi, bjet4_Phi);
		minDRbb6 = deltaR(bjet3_Eta, bjet4_Eta, bjet3_Phi, bjet4_Phi);
			
		minDRbb = minDRbb1;
		if(minDRbb2 < minDRbb) minDRbb = minDRbb2;
		if(minDRbb3 < minDRbb) minDRbb = minDRbb3;
		if(minDRbb4 < minDRbb) minDRbb = minDRbb4;
		if(minDRbb5 < minDRbb) minDRbb = minDRbb5;
		if(minDRbb6 < minDRbb) minDRbb = minDRbb6;
		
		Double_t minDRbj, minDRbj1, minDRbj2, minDRbj3, minDRbj4, minDRbj5, minDRbj6, minDRbj7, minDRbj8;
		
		minDRbj1 = deltaR(jet1_Eta, bjet1_Eta, jet1_Phi, bjet1_Phi);
		minDRbj2 = deltaR(jet1_Eta, bjet2_Eta, jet1_Phi, bjet2_Phi);
		minDRbj3 = deltaR(jet1_Eta, bjet3_Eta, jet1_Phi, bjet3_Phi);
		minDRbj4 = deltaR(jet1_Eta, bjet4_Eta, jet1_Phi, bjet4_Phi);
		minDRbj5 = deltaR(jet2_Eta, bjet1_Eta, jet2_Phi, bjet1_Phi);
		minDRbj6 = deltaR(jet2_Eta, bjet2_Eta, jet2_Phi, bjet2_Phi);
		minDRbj7 = deltaR(jet2_Eta, bjet3_Eta, jet2_Phi, bjet3_Phi);
		minDRbj8 = deltaR(jet2_Eta, bjet4_Eta, jet2_Phi, bjet4_Phi);
			
		minDRbj = minDRbj1;
		if(minDRbj2 < minDRbj) minDRbj = minDRbj2;
		if(minDRbj3 < minDRbj) minDRbj = minDRbj3;
		if(minDRbj4 < minDRbj) minDRbj = minDRbj4;
		if(minDRbj5 < minDRbj) minDRbj = minDRbj5;
		if(minDRbj6 < minDRbj) minDRbj = minDRbj6;
		if(minDRbj7 < minDRbj) minDRbj = minDRbj7;
		if(minDRbj8 < minDRbj) minDRbj = minDRbj8;
		
		Double_t DRjj = deltaR(jet1_Eta, jet2_Eta, jet1_Eta, jet2_Eta);
		Double_t DEjj = fabs(jet1_Eta - jet2_Eta);
		
		(SorB ? hDRjjS : hDRjjB)->Fill(DRjj, weight);
		(SorB ? hDRbjS : hDRbjB)->Fill(minDRbj, weight);
		(SorB ? hDRbbS : hDRbbB)->Fill(minDRbb, weight);
		(SorB ? hDEjjS : hDEjjB)->Fill(DEjj, weight);
		
		// Check kinematic requirements
		if(minDRbb < 0.75) continue;
		if(DEjj < 4.25) continue;
			
		tempKinCut += weight;
		tempErrorKinCut++;	
			
		// Set up temporary TLorentzVectors
		TLorentzVector test11, test12, test21, test22, test31, test32;
		test11 = vbJet1 + vbJet2;
		test12 = vbJet3 + vbJet4;
		test21 = vbJet1 + vbJet3;
		test22 = vbJet2 + vbJet4;
		test31 = vbJet1 + vbJet4;
		test32 = vbJet3 + vbJet2;
		
		Int_t nBJetPairsHiggsMass=0;
		
		if((test11.M() > 90) && (test11.M() < 135) && (test12.M() > 90) && (test12.M() < 135)) nBJetPairsHiggsMass++; 
		if((test21.M() > 90) && (test21.M() < 135) && (test22.M() > 90) && (test22.M() < 135)) nBJetPairsHiggsMass++;
		if((test31.M() > 90) && (test31.M() < 135) && (test32.M() > 90) && (test32.M() < 135)) nBJetPairsHiggsMass++;
		
		(SorB ? hMHHS : hMHHB)->Fill(vDiHiggs.M(), weight);
		(SorB ? hMjjS : hMjjB)->Fill(vDiJet.M(), weight);
		(SorB ? hMbbS : hMbbB)->Fill(nBJetPairsHiggsMass, weight);
		
		if(vDiHiggs.M() < 275) continue;
		if(nBJetPairsHiggsMass!=1) continue;
		if(vDiJet.M() < 400) continue;
		
		tempMassCut += weight;
		tempErrorMassCut++;
		
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
	
	histogram(hNLeptonsS, hNLeptonsB, c1, "Number of leptons", "Portion", "./Histograms/NLeptons_" + BackgroundSample + ".jpg");
	histogram(hNJetsPUS, hNJetsPUB, c1, "Number of non-pileup jets", "Portion", "./Histograms/NJetsPU_" + BackgroundSample + ".jpg");
	histogram(hNBJetsS, hNBJetsB, c1, "Number of bjets between light jets", "Portion", "./Histograms/NBJets_" + BackgroundSample + ".jpg");
	
	histogram(hDRjjS, hDRjjB, c1, "#Delta R between light jets", "Portion", "./Histograms/DRjj_" + BackgroundSample + ".jpg");
	histogram(hDRbjS, hDRbjB, c1, "Minimum #Delta R between bjets and light jets", "Portion", "./Histograms/DRbj_" + BackgroundSample + ".jpg");
	histogram(hDRbbS, hDRbbB, c1, "Minimum #Delta R between bjets", "Portion", "./Histograms/DRbb_" + BackgroundSample + ".jpg");
	histogram(hDEjjS, hDEjjB, c1, "#Delta #eta between light jets", "Portion", "./Histograms/DEjj_" + BackgroundSample + ".jpg");
	
	histogram(hMHHS, hMHHB, c1, "M_{HH}", "Portion", "./Histograms/MHH_" + BackgroundSample + ".jpg");
	histogram(hMjjS, hMjjB, c1, "M_{jj}", "Portion", "./Histograms/Mjj_" + BackgroundSample + ".jpg");
	histogram(hMbbS, hMbbB, c1, "Number of combinations with both M_{bb} in the range 90-135", "Portion", "./Histograms/Mbb_" + BackgroundSample + ".jpg");
	
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

	histogram(hNLeptonsS, c1, "Number of leptons", "Portion", "./Histograms/NLeptons_" + BackgroundSample + ".jpg");
	histogram(hNJetsPUS, c1, "Number of non-pileup jets", "Portion", "./Histograms/NJetsPU_" + BackgroundSample + ".jpg");
	histogram(hNBJetsS, c1, "Number of bjets between light jets", "Portion", "./Histograms/NBJets_" + BackgroundSample + ".jpg");
	
	histogram(hDRjjS, c1, "#Delta R between light jets", "Portion", "./Histograms/DRjj_" + BackgroundSample + ".jpg");
	histogram(hDRbjS, c1, "Minimum #Delta R between bjets and light jets", "Portion", "./Histograms/DRbj_" + BackgroundSample + ".jpg");
	histogram(hDRbbS, c1, "Minimum #Delta R between bjets", "Portion", "./Histograms/DRbb_" + BackgroundSample + ".jpg");
	histogram(hDEjjS, c1, "#Delta #eta between light jets", "Portion", "./Histograms/DEjj_" + BackgroundSample + ".jpg");
	
	histogram(hMHHS, c1, "M_{HH}", "Portion", "./Histograms/MHH_" + BackgroundSample + ".jpg");
	histogram(hMjjS, c1, "M_{jj}", "Portion", "./Histograms/Mjj_" + BackgroundSample + ".jpg");
	histogram(hMbbS, c1, "Number of combinations with both M_{bb} in the range 90-135", "Portion", "./Histograms/Mbb_" + BackgroundSample + ".jpg");
	
	// Print event yields
	cout << "Signal" << endl << endl;
	cout << "Total events: " << 2.01*3000*0.577*0.577 << endl;
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

/*
 * FUNCTION FOR dR calculation
 */
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
