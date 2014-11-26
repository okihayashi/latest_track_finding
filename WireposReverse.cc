#include <iostream>
#include <math.h>
#include "LayerInf140328.hh"
#include "WireposReverse.hh"

#define PI 3.14159265358979

void WireposReverse(double x, double y, int* layerID, double* theta){
    
    double radius = sqrt(x*x + y*y);
    //if((x > 0 && y > 0) || (x > 0 && y < 0) || (x < 0 && y > 0)){
    if(x > 0){
        *theta = atan(y/x);
    }
    if(x < 0){
	*theta = atan(y/x) + PI;
    }
    if(x == 0 && y != 0){
	*theta = acos(x/y)/asin(x/y);
    }
    if(x == 0 && y == 0){
	*theta = 0;
    }
    //*theta = asin(y/sqrt(x*x + y*y));

    double delta = 0.0001;
    for(int i=1;i<18;i++){
	LayerInf layerinf(i);
	double layer_radius = layerinf.GetR_0();
	if(fabs(radius - layer_radius) <= delta){
	    *layerID = i;
	    break;
	}
    }
}     
    
