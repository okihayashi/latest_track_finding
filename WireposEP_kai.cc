//-----------------------------------------------------------//                
//WireposEP.cc                                                                 
//                                                                             
//This function is used to get sense wires position (x,y,z)
//at end plete by using LayerInf140104.cc
//
//Argument: layer number of sense wire,
//	    cell number of sense wire,
//	    x-coordinate of position of sense wire,
//	    y-coordinate of position of sense wire,
//	    (z-coordinate of position of sense wire)  
//-----------------------------------------------------------//	               

#include <math.h>
#include "WireposEP_kai.hh"
#include "LayerInf140328.hh"

void WireposEP(int layer, double theta, double *x, double *y, double *z){

    double Radius = 0;
    double Radian = 0;
    int N;

    for(int i=0;i<18;i++){
	if(layer == i){
	    LayerInf layerinf(i);
	    Radius = layerinf.GetR_EP();
	    N = layerinf.GetNOfWire();
            *z = layerinf.Getz_EP();
	    double NOfSkip = static_cast<double>(layerinf.GetNOfSkip());
	    //Radian = theta + layerinf.GetInterval() * NOfSkip/2.;
	    Radian = theta - layerinf.GetInterval() * NOfSkip/2.; 
	    //for(int j=0;j<N;j++){
	    //    if(cell == j){
	    //        Radian = layerinf.GetInterval() * j;
	    //        break;
	    //    }
	    //}
	    //break;    
	}
    }
    *x = Radius * cos(Radian);
    *y = Radius * sin(Radian);
}
