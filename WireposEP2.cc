#include <math.h>
#include "WireposEP2.hh"
#include "LayerInf.hh"
#include <iostream>

using namespace std;
void WireposEP2(int layer, int cell, double *x, double *y, double *z){

    double Radius = 0;
    //double theta;
    double Radian = 0;
    int N;
    

    for(int i=0;i<20;i++){
	if(layer == i){
	    LayerInf layerinf(i); 
	    Radius = layerinf.GetR_EP();
	    N = layerinf.GetNOfWire();
	    *z = layerinf.Getz_EP() * (-1);
            //cout << "Skip :" << layerinf.GetNOfSkip() << endl; 

	    for(int j=0;j<N;j++){
		if(cell == j){
		    Radian = layerinf.GetInterval() * (j + layerinf.GetNOfSkip());
		    break;
		}
	    }
	}
    }

    *x = Radius * cos(Radian);
    *y = Radius * sin(Radian);

    
}
		    
