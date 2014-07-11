//-------------------------------------------------------------------
// Select bbtautau events from delphes
//
// execute with:
// root -l -q selectDelphes.C+\(\"_inputfile_name_\",_file_cross_section_\)
//
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

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

int puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar);

void selectDelphes(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToGGBB_14TeV/HHToGGBB_14TeV_0.root",
const Double_t xsec=1341.36923) {

	const TString inputFile = inputfile + ".txt";
	ifstream ifs(inputFile);

	assert(ifs.is_open());
	
	TString filename;
	
	cout << "Reading " << inputfile << endl;

	while(ifs >> filename){


	// read input input file
	TChain chain("Delphes");
	chain.Add("root://eoscms.cern.ch/" + (inputfile == "HHToBBBB_14TeV" ? "/store/user/arapyan/Delphes_phase2/VBFHHTobbbb_TuneZ2_14TeV-madgraph_june9/files/" : "/store/group/upgrade/delphes/ProdJun14/" + inputfile  + "/") + filename);
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// set up branches to read in from file
	TClonesArray *branchJet = treeReader->UseBranch("Jet");

	// set up loop variables
	Jet *jet;

	// set up storage variables
	Jet *bJet1=0, *bJet2=0, *bJet3=0, *bJet4=0, *Jet1=0, *Jet2=0;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;
	Double_t weight=3000000*xsec;
	UInt_t nEvents=0;

	const TString outfile = "/afs/cern.ch/work/a/ariostas/private/vbf-bbbb_temp/" + inputfile + "/" + filename;

	TFile *outFile = new TFile(outfile, "RECREATE");

	// tree to hold the number of events in the file before selection
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
	nEvents=numberOfEntries;
	sampTree->Fill();

	// tree to hold information about selected events
	TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("weight",		&weight,    "weight/D");  // number of jets
	outTree->Branch("bJet1",   		"Jet", &bJet1); // 4-vector for leading jet
	outTree->Branch("bJet2",   		"Jet", &bJet2); // 4-vector for second jet
	outTree->Branch("bJet3",   		"Jet", &bJet3); // 4-vector for second jet
	outTree->Branch("bJet4",   		"Jet", &bJet4); // 4-vector for second jet
	outTree->Branch("Jet1",   		"Jet", &Jet1); // 4-vector for second jet
	outTree->Branch("Jet2",   		"Jet", &Jet2); // 4-vector for second jet


	for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
		treeReader->ReadEntry(iEntry);

		// Reset index variables
		nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;
		
		// Select bjets and light jets
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
		
		// Check if there are four bjets and two light jets
		if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)) outTree->Fill();

	} // end event loop

	// save file
	outFile->Write();
	// close file
	outFile->Close();

	cout << "----SUMMARY----" << endl;
	cout << " input file " << filename << " selection done " << endl;

	}

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
