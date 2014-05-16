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
void analyzeB();
void saveResults();
Float_t deltaR( const Float_t, const Float_t, const Float_t, const Float_t);

Double_t Signal1=0, Signal2=0, Signal3=0, Signal4=0;

TH1D *hMJS = new TH1D("MJS", "MJS", 50, 0, 1000);



void diHiggs(){
	
	analyzeS();
	//analyzeB();
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
	Signal1 = numberOfEntries;

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
				
				Signal2++;
				
				if((Jet1->PT >30) && (Jet2->PT >30) && (bJet1->PT >30) && (bJet2->PT >30) && (bJet3->PT >30) && (bJet4->PT >30) &&
					(fabs(Jet1->Eta)<5) && (fabs(Jet2->Eta)<5) && (fabs(bJet1->Eta)<5) && (fabs(bJet2->Eta)<5) && (fabs(bJet3->Eta)<5) && (fabs(bJet4->Eta)<5)){
					
					Signal3++;
					
					vJet1.SetPtEtaPhiM(Jet1->PT, Jet1->Eta, Jet1->Phi, Jet1->Mass);
					vJet2.SetPtEtaPhiM(Jet2->PT, Jet2->Eta, Jet2->Phi, Jet2->Mass);
					
					vJetSys = vJet1 + vJet2;
					
					hMJS->Fill(vJetSys.M());
					
					if((vJetSys.M() > 500) && (fabs(vJetSys.Eta()) > 2.5)){
						
						Signal4++;
						
					}
					
				}
				
			}
		
	} // end event loop

	//preSignalEvents = preSignalEvents*3000*40*0.213*0.213/numberOfEntries;
	//signalEvents = signalEvents*3000*40*0.213*0.213/numberOfEntries;
	
	ifs.close();
}


void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 900, 600);
	
	gStyle->SetOptStat(kFALSE);

	histogram(hMJS, c1, "Mass diJet system", "Count", "histogramL1Pt.jpg");
	
	cout << "Total events: " << 1.8*3000*0.577*0.577 << "\nEvents with all jets: " << 1.8*3000*0.577*0.577*Signal2/Signal1 << endl;
	cout << "Events with selected jets: " << 1.8*3000*0.577*0.577*Signal3/Signal1 << "\nEvents after dijet cut: " << 1.8*3000*0.577*0.577*Signal4/Signal1 << endl;
	
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
