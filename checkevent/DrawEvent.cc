#include <iostream> 
#include <vector>
#include <TTree.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#define PI 3.14159265358979

using namespace std;

class LayerInf {
    private:
	int NOfWire;
	double Interval;
	double R_0;
	double NOfSkip;

    public:
	LayerInf(int l);
	
	double GetNOfWire()		{ return NOfWire;  }
	double GetInterval()		{ return Interval; }
	double GetR_0()			{ return R_0;	   }
	double GetNOfSkip()		{ return NOfSkip;  }
};

LayerInf::LayerInf(int l){
    double r0[20] = {51.1673,52.7600,54.3671,55.9738,57.5801,59.1860,60.7917,62.3970,64.0021,65.6069,                    
        	     67.2115,68.8159,70.4201,72.0241,73.6280,75.1709,76.7756,78.3802,79.9846,81.5889}; //[cm]

    int numberofwire[20] = {396/2,396/2,408/2,420/2,432/2,444/2,456/2,468/2,480/2,492/2,  
        		    504/2,516/2,528/2,540/2,552/2,564/2,576/2,588/2,600/2,612/2}; 

    double interval[20];                    
    for(int i=0;i<20;i++){                  
        interval[i] = 2*PI/numberofwire[i]; 
    }

    int NOfskippinghole[20] = {6,-6,6,-6,6,-6,6,-6,6,-6,6,-6,6,-6,6,-7,7,-7,7,-7};      

    NOfWire = numberofwire[l];    
    Interval = interval[l];       
    R_0 = r0[l];                  
    NOfSkip = NOfskippinghole[l]; 
}

void Wirepos0(int layer,int cell,double* x, double* y){

   double Radius;                                                                          
   double Radian;                                                                          
   int N;                                                                                  
                                                                                           
   for(int i=0;i<20;i++){                                                                  
    if(layer == i){                                                                     
   	LayerInf layerinf(i);                                                           
   	Radius = layerinf.GetR_0();                                                     
           N = layerinf.GetNOfWire();                                                      
   	for(int j=0;j<N;j++){                                                           
   	    if(cell == j){                                                              
   		Radian = layerinf.GetInterval()*(j + double(layerinf.GetNOfSkip()/2.)); 
   		break;                                                                  
   	    }                                                                           
   	}                                                                               
    }                                                                                   
   }                                                                                       
   *x = Radius * cos(Radian);                                                              
   *y = Radius * sin(Radian);                                                              
}                                                                                           

void DrawEvent(){

    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    TFile* file = new TFile("../../Data/signal.140905M02.noise-3.root");
    TTree* t = (TTree*)file->Get("tree");

    TCanvas* c1 = new TCanvas("c1","c1",10,10,800,800);

    int CDCcell_nHits;
    vector<double>* CDCcell_edep = 0;
    vector<double>* CDCcell_layerID = 0;
    vector<double>* CDCcell_cellID = 0;
    vector<double>* CDCcell_pz = 0;
    vector<int>* CDCcell_hittype = 0;

    t->SetBranchAddress("CdcCell_nHits",&CDCcell_nHits);
    t->SetBranchAddress("CdcCell_edep",&CDCcell_edep);
    t->SetBranchAddress("CdcCell_layerID",&CDCcell_layerID);
    t->SetBranchAddress("CdcCell_cellID",&CDCcell_cellID);
    t->SetBranchAddress("CdcCell_pz",&CDCcell_pz);
    t->SetBranchAddress("CdcCell_hittype",&CDCcell_hittype);
    
    //for(int a=0;a<t->GetEntries();i++){
    for(int a=0;a<1;a++){
	t->GetEntry(1);
	
	int count_signal = 0;
	for(int i=0;i<CDCcell_nHits;i++){
	    if(CDCcell_hittype->at(i) == 0){
		count_signal++;
	    }
	}
	if(count_signal>14){
	    vector<double> CDCcell_pz_sig;
	    vector<double> CDCcell_x;
	    vector<double> CDCcell_y;
	    
	    for(int i=0;i<CDCcell_nHits;i++){
		if(CDCcell_edep->at(i)*1e6 > 7.5){
		    continue;
		}else{
		    double x,y;
		    Wirepos0(CDCcell_layerID->at(i),CDCcell_cellID->at(i),&x,&y);
		    CDCcell_x.push_back(x);
		    CDCcell_y.push_back(y);
                    if(CDCcell_hittype->at(i) == 0){
			CDCcell_pz_sig.push_back(CDCcell_pz->at(i));
		    }
		}
	    }
	    break;
	}
    }
    TGraph* g1 = new TGraph(CDCcell_x.size(),&CDCcell_x[0],&CDCcell_y[0]);
    g1->Draw("ap");
    g1->SetMarkerStyle(4);
    g1->SetMarkerSize(0.2);
    g1->SetMarkerColor(4);
		
    c1->Update();
}    

	    

	    

