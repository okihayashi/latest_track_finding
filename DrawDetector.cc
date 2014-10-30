//-----------------------------------------------------------//         
//DrawDetector.cc                                                        
//                                                                      
//This function is used to draw all sense wires in 
//event display                 
//-----------------------------------------------------------//	        

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "Wirepos0.hh"
#include "LayerInf140328.hh"
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>

using namespace std;
void DrawDetector(){

    double x_ep, y_ep;
    double x[18][300];
    double y[18][300];
    TGraph *g[18];
    for(int i=0;i<18;i++){
	LayerInf layerinf(i+1);
        int N = layerinf.GetNOfWire();
	for(int j=0;j<N;j++){
	    x_ep = 0;
            y_ep = 0;
	    Wirepos0(i+1, j, &x_ep, &y_ep);
            x[i][j] = x_ep;
	    y[i][j] = y_ep;
	}
	g[i] = new TGraph(N, x[i], y[i]);
	if(i == 0){
	    g[i]->Draw("p,same");
            g[i]->SetMarkerStyle(5);  
            g[i]->SetMarkerSize(0.2); 
            g[i]->SetMarkerColor(17); 
	    //g[i]->GetXaxis()->SetLimits(-90,90);
	    //g[i]->SetMinimum(-90);
	    //g[i]->SetMaximum(90);
	}else{
	    g[i]->Draw("p,same");
	    g[i]->SetMarkerStyle(5);   
            g[i]->SetMarkerSize(0.2);  
	    g[i]->SetMarkerColor(17);  
	}
    }
}

/*int main(int argc, char** argv){
    
    TApplication app("app",&argc,argv); 
    TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800); 
    c1->cd();
    DrawDetector();

    app.Run();
}
*/
