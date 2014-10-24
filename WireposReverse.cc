#include <iostream>
#include <math.h>
#include "LayerInf140328.hh"
#include "WireposReverse.hh"

void WireposReverse(double x, double y, int* layerID, double* theta){
    
    double radius = sqrt(x*x + y*y);
    *theta = acos(x/sqrt(x*x + y*y));
    

    for(int i=0;i<20;i++){
	LayerInf layerinf(i);
	double layer_radius = layerinf.GetR_0();
	if(layer_radius == radius){
	    *layerID = i;
	    break;
	}
    }
}     
    
