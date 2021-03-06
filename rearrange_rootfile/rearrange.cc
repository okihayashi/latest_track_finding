#include <stdio.h>
#include <stdlib.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <vector>
#include <iostream>

using namespace std;
int main(int argc, char** argv){

    TFile* file = new TFile("/Users/hayashi_oki/Workspace/Track-finding/Helix/NN/Data/signal.140905M02.noise-2.root");
    TTree* t = (TTree*)file->Get("tree");

    int CdcCell_nHits = 0;
    vector<double>* CdcCell_x = 0;
    vector<double>* CdcCell_y = 0;
    vector<double>* CdcCell_z = 0;
    vector<double>* CdcCell_px = 0;
    vector<double>* CdcCell_py = 0;
    vector<double>* CdcCell_pz = 0;
    vector<double>* CdcCell_cellID = 0;
    vector<double>* CdcCell_layerID = 0;
    vector<double>* CdcCell_edep = 0;
    vector<int>* CdcCell_hittype = 0;

    t->SetBranchAddress("CdcCell_nHits",&CdcCell_nHits);
    t->SetBranchAddress("CdcCell_px",&CdcCell_px);
    t->SetBranchAddress("CdcCell_x",&CdcCell_x);
    t->SetBranchAddress("CdcCell_y",&CdcCell_y);
    t->SetBranchAddress("CdcCell_z",&CdcCell_z);
    t->SetBranchAddress("CdcCell_py",&CdcCell_py);
    t->SetBranchAddress("CdcCell_pz",&CdcCell_pz);
    t->SetBranchAddress("CdcCell_cellID",&CdcCell_cellID);
    t->SetBranchAddress("CdcCell_layerID",&CdcCell_layerID);
    t->SetBranchAddress("CdcCell_edep",&CdcCell_edep);
    t->SetBranchAddress("CdcCell_hittype",&CdcCell_hittype);

    TFile* newfile = TFile::Open("/Users/hayashi_oki/Workspace/Track-finding/Helix/NN/Data/signal.140905M02.noise-3.root","RECREATE");
    TTree* t2 = new TTree("tree","tree");
    
    int CDCcell_eventID = 0;
    int CDCcell_nHits = 0;
    vector<double> CDCcell_x(0);
    vector<double> CDCcell_y(0);
    vector<double> CDCcell_z(0);
    vector<double> CDCcell_px(0);      
    vector<double> CDCcell_py(0);      
    vector<double> CDCcell_pz(0);      
    vector<double> CDCcell_cellID(0);  
    vector<double> CDCcell_layerID(0); 
    vector<double> CDCcell_edep(0);    
    vector<int> CDCcell_hittype(0);   
    
    t2->Branch("CdcCell_eventID",&CDCcell_eventID);
    t2->Branch("CdcCell_nHits",&CDCcell_nHits);
    t2->Branch("CdcCell_x",&CDCcell_x);
    t2->Branch("CdcCell_y",&CDCcell_y);
    t2->Branch("CdcCell_z",&CDCcell_z);
    t2->Branch("CdcCell_px",&CDCcell_px);
    t2->Branch("CdcCell_py",&CDCcell_py);
    t2->Branch("CdcCell_pz",&CDCcell_pz);
    t2->Branch("CdcCell_cellID",&CDCcell_cellID);
    t2->Branch("CdcCell_layerID",&CDCcell_layerID);
    t2->Branch("CdcCell_edep",&CDCcell_edep);
    t2->Branch("CdcCell_hittype",&CDCcell_hittype);

    int flag;

    for(int i=0;i<t->GetEntries();i++){
	cout << "i = " << i << endl;
	t->GetEntry(i);
        CDCcell_nHits = CdcCell_nHits;

        for(int a=0;a<CDCcell_nHits;a++){
	    CDCcell_x.push_back(CdcCell_x->at(a));
	    CDCcell_y.push_back(CdcCell_y->at(a)); 
            CDCcell_z.push_back(CdcCell_z->at(a)); 
	    CDCcell_px.push_back(CdcCell_px->at(a));
            CDCcell_py.push_back(CdcCell_py->at(a)); 
            CDCcell_pz.push_back(CdcCell_pz->at(a));
            CDCcell_cellID.push_back(CdcCell_cellID->at(a));
            CDCcell_layerID.push_back(CdcCell_layerID->at(a));
	    CDCcell_edep.push_back(CdcCell_edep->at(a));
	    CDCcell_hittype.push_back(CdcCell_hittype->at(a));
	}

        double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy14, dummy15, dummy16;
	do{
	    flag = 0;
	    for(int j=0;j<CDCcell_nHits-1;j++){
		if(CDCcell_layerID[j] > CDCcell_layerID[j+1]){
		    //cout << "j = " << j << endl;
		    flag = 1;
                    
		    dummy1 = CDCcell_px[j];
		    CDCcell_px[j] = CDCcell_px[j+1];
		    CDCcell_px[j+1] = dummy1;

                    dummy2 = CDCcell_py[j];              
                    CDCcell_py[j] = CDCcell_py[j+1]; 
                    CDCcell_py[j+1] = dummy2;

                    dummy3 = CDCcell_pz[j];               
                    CDCcell_pz[j] = CDCcell_pz[j+1];  
                    CDCcell_pz[j+1] = dummy3;             

                    dummy4 = CDCcell_cellID[j];               
                    CDCcell_cellID[j] = CDCcell_cellID[j+1];  
                    CDCcell_cellID[j+1] = dummy4;             

                    dummy5 = CDCcell_layerID[j];                    
                    CDCcell_layerID[j] = CDCcell_layerID[j+1];   
                    CDCcell_layerID[j+1] = dummy5;

                    dummy6 = CDCcell_edep[j];                    
                    CDCcell_edep[j] = CDCcell_edep[j+1];   
                    CDCcell_edep[j+1] = dummy6;

                    dummy7 = CDCcell_hittype[j];                    
                    CDCcell_hittype[j] = CDCcell_hittype[j+1];   
                    CDCcell_hittype[j+1] = dummy7;

		    dummy14 = CDCcell_x[j];                     
                    CDCcell_x[j] = CDCcell_x[j+1];       
                    CDCcell_x[j+1] = dummy14;                   

                    dummy15 = CDCcell_y[j];                      
                    CDCcell_y[j] = CDCcell_y[j+1];               
                    CDCcell_y[j+1] = dummy15;                    

                    dummy16 = CDCcell_z[j];                      
                    CDCcell_z[j] = CDCcell_z[j+1];               
                    CDCcell_z[j+1] = dummy16;                    
		}
	    }
	}while(flag == 1);

        //for(int k=0;k<CDCcell_nHits;k++){
        //    cout << "CDCcell_edep[" << k << "] = " << CDCcell_edep[k] << endl;
	//}

	int layerID;
        int N[18] = {};
	double dummy8, dummy9, dummy10, dummy11, dummy12, dummy13, dummy17, dummy18, dummy19;

	for(int i=0;i<CDCcell_nHits;i++){
	    layerID = CDCcell_layerID[i];
	    N[layerID]++;
	}

        do{
            int n = -1;
            int start = 0;
	    flag = 0;
	    for(int i=0;i<18;i++){
		n += N[i];
		for(int j=0;j<CDCcell_nHits;j++){
		    if(CDCcell_layerID[j] == i+1){
			for(int k=start;k<n;k++){
			    if(CDCcell_cellID[k] > CDCcell_cellID[k+1]){
				flag = 1;
				
                                dummy8 = CDCcell_px[k];                           
                                CDCcell_px[k] = CDCcell_px[k+1];              
                                CDCcell_px[k+1] = dummy8;                         
                                                                                      
                                dummy9 = CDCcell_py[k];                           
                                CDCcell_py[k] = CDCcell_py[k+1];              
                                CDCcell_py[k+1] = dummy9;                         
                                                                                      
                                dummy10 = CDCcell_pz[k];                           
                                CDCcell_pz[k] = CDCcell_pz[k+1];              
                                CDCcell_pz[k+1] = dummy10;                         
                                                                                      
                                dummy11 = CDCcell_cellID[k];                       
                                CDCcell_cellID[k] = CDCcell_cellID[k+1];      
                                CDCcell_cellID[k+1] = dummy11;                     
                                                                                      
                                dummy12 = CDCcell_edep[k];                         
                                CDCcell_edep[k] = CDCcell_edep[k+1];          
                                CDCcell_edep[k+1] = dummy12;                       
                                                                                      
                                dummy13 = CDCcell_hittype[k];                      
                                CDCcell_hittype[k] = CDCcell_hittype[k+1];    
                                CDCcell_hittype[k+1] = dummy13;
			    
			        dummy17 = CDCcell_x[k];       
			        CDCcell_x[k] = CDCcell_x[k+1];
			        CDCcell_x[k+1] = dummy17;     
			    
			        dummy18 = CDCcell_y[k];       
			        CDCcell_y[k] = CDCcell_y[k+1];
			        CDCcell_y[k+1] = dummy18;     
			    
			        dummy19 = CDCcell_z[k];       
                                CDCcell_z[k] = CDCcell_z[k+1];
                                CDCcell_z[k+1] = dummy19;     
			    }
			}
		    }
		}
		start += N[i];
	    }
	}while(flag == 1);

	t2->Fill();
	CDCcell_eventID++;
	CDCcell_x.clear();
	CDCcell_y.clear();
	CDCcell_z.clear();
        CDCcell_px.clear();      
        CDCcell_py.clear();      
        CDCcell_pz.clear();      
        CDCcell_cellID.clear();  
        CDCcell_layerID.clear(); 
        CDCcell_edep.clear();    
        CDCcell_hittype.clear(); 
    }
    delete file;
    newfile->cd();
    t2->Write();
    delete newfile;
}

