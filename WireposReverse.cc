#include <iostream>
#include <math.h>
#include "LayerInf140328.hh"
#include "WireposReverse.hh"

void WireposReverse(double x, double y, int* layerID, double* theta){
    
    double radius = sqrt(x*x + y*y);
    
    if(x != 0){
	*theta = atan(y/x);
    }
    if(x == 0 && y != 0){
	*theta = acos(x/y)/asin(x/y);
    }
    if(x == 0 && y == 0){
	*theta = 0;
    }
    
    //*theta = asin(y/sqrt(x*x + y*y));

    double delta = 0.0001;
    for(int i=0;i<20;i++){
	LayerInf layerinf(i);
	double layer_radius = layerinf.GetR_0();
	if(fabs(radius - layer_radius) <= delta){
	    *layerID = i;
	    break;
	}
    }
}     
    
