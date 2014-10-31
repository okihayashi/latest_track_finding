#include <math.h>
#include "WireposEP2_kai.hh"
#include "LayerInf140328.hh"

void WireposEP2(int layer, double theta, double *x, double *y, double *z){

    double Radius = 0;
    double Radian = 0;
    int N;
    for(int i=0;i<20;i++){
	if(layer == i){
	    LayerInf layerinf(i); 
	    Radius = layerinf.GetR_EP();
	    N = layerinf.GetNOfWire();
	    *z = layerinf.Getz_EP();
            double NOfSkip = static_cast<double>(layerinf.GetNOfSkip());
	    Radian = theta + layerinf.GetInterval() * NOfSkip/2.;
	    //for(int j=0;j<N;j++){
	    //    if(cell == j){
	    //        Radian = layerinf.GetInterval() * (j + layerinf.GetNOfSkip());
	    //        break;
	    //    }
	    //}
	}
    }
    *x = Radius * cos(Radian);
    *y = Radius * sin(Radian);
}
		    
