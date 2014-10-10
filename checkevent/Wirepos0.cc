#include <math.h>
#include "Wirepos0.hh"
#include "LayerInf140328.hh"

void Wirepos0(int layer, int cell, double *x, double *y){

    double Radius;
    double Radian;
    int N;
    
    //--- Return Radius and Theta(Radian) at Z=0
    for(int i=0;i<20;i++){
	if(layer == i){
	    LayerInf layerinf(i);
	    Radius = layerinf.GetR_0();
            N = layerinf.GetNOfWire();
	    for(int j=0;j<N;j++){
		if(cell == j){
		    Radian = layerinf.GetInterval()*j;
		    break;
		}
	    }
	}
    }
    //--- Calculate (x, y) at Z=0
    *x = Radius * cos(Radian);
    *y = Radius * sin(Radian);
}








