#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TAxis.h>
#include "Wirepos0.hh"

using namespace std;
int main(int argc, char** argv){

    TApplication app("app",&argc,argv);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    TFile* file = new TFile("../../Data/signal.140905M02.noise-3.root");
    TTree* t = (TTree*)file->Get("tree");

    TCanvas* c1 = new TCanvas("c1","c1",10,10,800,800);

    int CDCcell_nHits;
    vector<double>* CDCcell_edep = 0;
    vector<double>* CDCcell_layerID = 0;
    vector<double>* CDCcell_cellID = 0;
    vector<int>*    CDCcell_hittype = 0;

    t->SetBranchAddress("CdcCell_nHits",&CDCcell_nHits); 
    t->SetBranchAddress("CdcCell_edep",&CDCcell_edep);
    t->SetBranchAddress("CdcCell_layerID",&CDCcell_layerID);
    t->SetBranchAddress("CdcCell_cellID",&CDCcell_cellID);
    t->SetBranchAddress("CdcCell_hittype",&CDCcell_hittype);

    //for(int k=0;k<1000;k++){
    for(int k=0;k<1;k++){
        t->GetEntry(1);

        vector<double> CDCcell_x;
        vector<double> CDCcell_y;
        vector<double> CDCcell_signalhits_x;
        vector<double> CDCcell_signalhits_y;
	vector<double> CDCcell_bghits_x;
	vector<double> CDCcell_bghits_y;
	vector<double> CDCcell_hittype_cut;
	vector<double> CDCcell_layerID_cut;
	vector<double> CDCcell_cellID_cut;

	for(int i=0;i<CDCcell_nHits;i++){
	    if(CDCcell_edep->at(i)*1e6 > 7.5){   //keV
	        continue;
	    }else{
	        double x0;
		double y0;
		Wirepos0(CDCcell_layerID->at(i),CDCcell_cellID->at(i),&x0,&y0);
		CDCcell_x.push_back(x0);
	        CDCcell_y.push_back(y0);
		CDCcell_hittype_cut.push_back(CDCcell_hittype->at(i));
		CDCcell_layerID_cut.push_back(CDCcell_layerID->at(i));
		CDCcell_cellID_cut.push_back(CDCcell_cellID->at(i));
	    }
	}
        for(int i=0;i<CDCcell_x.size();i++){
	    if(CDCcell_hittype_cut.at(i) == 0){
		double x_signal;
		double y_signal;
                Wirepos0(CDCcell_layerID_cut.at(i),CDCcell_cellID_cut.at(i),&x_signal,&y_signal);
		CDCcell_signalhits_x.push_back(x_signal);
		CDCcell_signalhits_y.push_back(y_signal);
	    }else{
		double x_bg;
		double y_bg;
		Wirepos0(CDCcell_layerID_cut.at(i),CDCcell_cellID_cut.at(i),&x_bg,&y_bg);
		CDCcell_bghits_x.push_back(x_bg);
		CDCcell_bghits_y.push_back(y_bg);
	    }
	}

        cout << "CDCcell_signalhits_x.size() = " << CDCcell_signalhits_x.size() << endl;

	TGraph* g1 = new TGraph(CDCcell_x.size(),&CDCcell_x[0],&CDCcell_y[0]);
	c1->cd();
	g1->Draw("AP");
	g1->SetMarkerStyle(4);
	g1->SetMarkerSize(0.3);
	g1->SetMarkerColor(11);
	g1->SetMinimum(-90);
	g1->SetMaximum(90);
	g1->GetXaxis()->SetLimits(-90,90);

        TGraph* gsig = new TGraph(CDCcell_signalhits_x.size(),&CDCcell_signalhits_x[0],&CDCcell_signalhits_y[0]); 
        gsig->Draw("P,same");                                                        
        gsig->SetMarkerStyle(4);                                                 
        gsig->SetMarkerSize(0.3);                                                
        gsig->SetMarkerColor(2);

        TGraph* gbg = new TGraph(CDCcell_bghits_x.size(),&CDCcell_bghits_x[0],&CDCcell_bghits_y[0]);  
        gbg->Draw("P,same");                                                                                      
        gbg->SetMarkerStyle(4);                                                                                   
        gbg->SetMarkerSize(0.3);                                                                                  
        gbg->SetMarkerColor(4);                                                                                   


    }
    c1->Update();
    app.Run();
}
