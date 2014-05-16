#if !defined(__CINT__) || defined(__MAKECINT__)
// include statements for all needed dependencies
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

// include statements for Delphes
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void histogram(TH1D*, TH1D*, TCanvas, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void analyzeS();
void analyzeB(TString, Double_t);
void saveResults();
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);

Double_t totalSignal=0, jetsSignal=0, selJetsSignal=0, dijetCutSignal=0, injetCutSignal=0;
Double_t finalBG1=0, finalBG2=0, finalBG3=0, finalBG4=0, finalBG5=0;

TH1D *hMJSS = new TH1D("MJSB", "MJSB", 50, 0, 1000);
TH1D *hMJSB = new TH1D("MJSB", "MJSB", 50, 0, 1000);
TH1D *hInBJetsS = new TH1D("InJetsS", "InJetsS", 6, -0.5, 5.5);
TH1D *hInBJetsB = new TH1D("InJetsB", "InJetsB", 6, -0.5, 5.5);



void diHiggs(){
	
	analyzeS();
	analyzeB("tt-4p-0-600-v1510_14TEV", 530.89358);
	analyzeB("tt-4p-600-1100-v1510_14TEV", 42.55351);
	analyzeB("tt-4p-1100-1700-v1510_14TEV", 4.48209);
	analyzeB("tt-4p-1700-2500-v1510_14TEV", 0.52795);
	analyzeB("tt-4p-2500-100000-v1510_14TEV", 0.05449);
	saveResults();
}



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

	// set up loop variables
	Jet *jet;

	// set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4, *Jet1, *Jet2;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;

	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vJetSys;

	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");

	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // event loop
		treeReader->ReadEntry(iEntry);

		// reset counter variable
		nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;

			
			for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // jet loop
				jet = (Jet*) branchJet->At(iJet);
				
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
			
			
			
			
			if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)){
				
				jetsSignal++;
				
				if((Jet1->PT >30) && (Jet2->PT >30) && (bJet1->PT >30) && (bJet2->PT >30) && (bJet3->PT >30) && (bJet4->PT >30) &&
					(fabs(Jet1->Eta)<5) && (fabs(Jet2->Eta)<5) && (fabs(bJet1->Eta)<5) && (fabs(bJet2->Eta)<5) && (fabs(bJet3->Eta)<5) && (fabs(bJet4->Eta)<5)){
					
					selJetsSignal++;
					
					vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
					vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
					
					vJetSys = vJet1 + vJet2;
					
					hMJSS->Fill(vJetSys.M());
					
					if((vJetSys.M() > 500) && (fabs(Jet1->Eta - Jet2->Eta) > 2.5)){
						
						dijetCutSignal++;
						
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
						
						hInBJetsS->Fill(inBJets);
						
						if(inBJets==4)injetCutSignal++;
						
					}
					
				}
				
			}
		
	} // end event loop

	//preSignalEvents = preSignalEvents*3000*40*0.213*0.213/numberOfEntries;
	//signalEvents = signalEvents*3000*40*0.213*0.213/numberOfEntries;
	
	ifs.close();
}

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

	// set up loop variables
	Jet *jet;

	// set up storage variables
	Jet *bJet1, *bJet2, *bJet3, *bJet4, *Jet1, *Jet2;
	Int_t nbJet1, nbJet2, nbJet3, nbJet4, nJet1, nJet2;

	TLorentzVector vbJet1, vbJet2, vbJet3, vbJet4, vJet1, vJet2, vJetSys;

	
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchEvent = treeReader->UseBranch("Event");

	for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // event loop
		treeReader->ReadEntry(iEntry);

		// reset counter variable
		nbJet1=nbJet2=nbJet3=nbJet4=nJet1=nJet2=-1;

			
			for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // jet loop
				jet = (Jet*) branchJet->At(iJet);
				
				if(jet->BTag){
					if(nbJet1==-1){
						nbJet1=iJet;
						bJet1= (Jet*) branchJet->At(nbJet1);
					}
					else if(jet->PT > bJet1->PT){
						
						nbJet4=nbJet3;
						bJet4= (Jet*) branchJet->At(nbJet4);
						
						nbJet3=nbJet2;
						bJet3= (Jet*) branchJet->At(nbJet3);
						
						nbJet2=nbJet1;
						bJet2= (Jet*) branchJet->At(nbJet2);
						
						nbJet1=iJet;
						bJet1= (Jet*) branchJet->At(nbJet1);
						
					}
					else if(nbJet2==-1){
						
						nbJet2=iJet;
						bJet2= (Jet*) branchJet->At(nbJet2);
						
					}
					else if(jet->PT > bJet2->PT){
						
						nbJet4=nbJet3;
						bJet4= (Jet*) branchJet->At(nbJet4);
						
						nbJet3=nbJet2;
						bJet3= (Jet*) branchJet->At(nbJet3);
						
						nbJet2=iJet;
						bJet2= (Jet*) branchJet->At(nbJet2);
						
					}
					else if(nbJet3==-1){
						
						nbJet3=iJet;
						bJet3= (Jet*) branchJet->At(nbJet3);
						
					}
					else if(jet->PT > bJet3->PT){
						
						nbJet4=nbJet3;
						bJet4= (Jet*) branchJet->At(nbJet4);
						
						nbJet3=nbJet2;
						bJet3= (Jet*) branchJet->At(nbJet3);
						
					}
					else if(nbJet4==-1){
						
						nbJet4=iJet;
						bJet4= (Jet*) branchJet->At(nbJet4);
						
					}
					else if(jet->PT > bJet4->PT){
						
						nbJet4=iJet;
						bJet4= (Jet*) branchJet->At(nbJet4);
						
					}
				}
				else{
					if(nJet1==-1){
						nJet1=iJet;
						Jet1= (Jet*) branchJet->At(nJet1);
					}
					else if((jet->PT > Jet1->PT)){
						nJet2=nJet1;
						Jet2= (Jet*) branchJet->At(nJet2);
						nJet1=iJet;
						Jet1= (Jet*) branchJet->At(nJet1);
						}
					else if(nJet2==-1){
						nJet2=iJet;
						Jet2= (Jet*) branchJet->At(nJet2);
					}
					else if(jet->PT > Jet2->PT){
						nJet2=iJet;
						Jet2= (Jet*) branchJet->At(nJet2);
					}
					
				}
			}
			
			
			
			
			if((nJet1!=-1) && (nJet2!=-1) && (nbJet1!=-1) && (nbJet2!=-1) && (nbJet3!=-1) && (nbJet4!=-1)){
				
				LHEFEvent *event = (LHEFEvent*) branchEvent->At(0);
				weight = event->Weight;
				
				Background2+=weight;
				
				if((Jet1->PT >30) && (Jet2->PT >30) && (bJet1->PT >30) && (bJet2->PT >30) && (bJet3->PT >30) && (bJet4->PT >30) &&
					(fabs(Jet1->Eta)<5) && (fabs(Jet2->Eta)<5) && (fabs(bJet1->Eta)<5) && (fabs(bJet2->Eta)<5) && (fabs(bJet3->Eta)<5) && (fabs(bJet4->Eta)<5)){
					
					Background3+=weight;
					
					vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
					vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
					
					vJetSys = vJet1 + vJet2;
					
					hMJSB->Fill(vJetSys.M());
					
					if((vJetSys.M() > 500) && (fabs(Jet1->Eta - Jet2->Eta) > 2.5)){
						
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
						
						hInBJetsB->Fill(inBJets);
						
						if(inBJets==4)Background5+=weight;
						
					}
					
				}
				
			}
		
	} // end event loop
	
	finalBG1 += crossSection*3000000*Background1/numberOfEntries;
	finalBG2 += crossSection*3000000*Background2/numberOfEntries;
	finalBG3 += crossSection*3000000*Background3/numberOfEntries;
	finalBG4 += crossSection*3000000*Background4/numberOfEntries;
	finalBG5 += crossSection*3000000*Background5/numberOfEntries;

	//preSignalEvents = preSignalEvents*3000*40*0.213*0.213/numberOfEntries;
	//signalEvents = signalEvents*3000*40*0.213*0.213/numberOfEntries;
	
	ifs.close();
}


void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 900, 600);
	
	gStyle->SetOptStat(kFALSE);

	histogram(hMJSS, c1, "Mass diJet system", "Count", "histogramMJSS.jpg");
	histogram(hMJSB, c1, "Mass diJet system", "Count", "histogramMJSB.jpg");
	histogram(hInBJetsS, c1, "Number of bjets between light jets", "Count", "histogramInBJetsS.jpg");
	histogram(hInBJetsB, c1, "Number of bjets between light jets", "Count", "histogramInBJetsB.jpg");
	
	cout << "Total events: " << 1.8*3000*0.577*0.577 << endl;
	cout << "Events with all jets: " << 1.8*3000*0.577*0.577*jetsSignal/totalSignal << endl;
	cout << "Events with selected jets: " << 1.8*3000*0.577*0.577*selJetsSignal/totalSignal << endl;
	cout << "Events after dijet cut: " << 1.8*3000*0.577*0.577*dijetCutSignal/totalSignal << endl;
	cout << "Events after injet cut: " << 1.8*3000*0.577*0.577*injetCutSignal/totalSignal << endl;
	
	cout << "\nBackground2: " << finalBG2 << endl;
	cout << "Background3: " << finalBG3 << endl;
	cout << "Background4: " << finalBG4 << endl;
	cout << "Background4: " << finalBG5 << endl;
	
	/*ofstream outputFile;
	outputFile.open ("/afs/cern.ch/user/a/ariostas/files9/histograms/logfile.txt");
	outputFile << "Process completed successfully" << endl;
	outputFile.close();*/

}


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

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

	const Float_t pi = 3.14159265358979;

	Float_t etaDiff = (eta1-eta2);
	Float_t phiDiff = fabs(phi1-phi2);
	while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

	Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

	return TMath::Sqrt(deltaRSquared);

}
