#include <math.h>
#include "Wirepos0.hh"
#include "LayerInf140328.hh"

void Wirepos0(int layer, int cell, double *x, double *y){

    double Radius;
    double Radian;
    int N;
    
    //--- Return Radius and Theta(Radian) at Z=0
    for(int i=0;i<18;i++){
	if(layer == i){
            double Radian_Offset = 0;
	    LayerInf layerinf(i);
            if(i==0 || i==1 || i==3 || i==6 || i==8 || i==10 || i==12 || i==16){
                Radian_Offset = layerinf.GetInterval()/2.;
            }
	    Radius = layerinf.GetR_0();
            N = layerinf.GetNOfWire();
	    double NOfSkip = static_cast<double>(layerinf.GetNOfSkip());
	    for(int j=0;j<N;j++){
		if(cell == j){
		    Radian = Radian_Offset + layerinf.GetInterval()*(j + NOfSkip/2.);
		    break;
		}
	    }
        }
    }
    //--- Calculate (x, y) at Z=0
    *x = Radius * cos(Radian);
    *y = Radius * sin(Radian);
}








