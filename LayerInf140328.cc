//-----------------------------------------------------------//            
//LayerInf.cc                                                              
//                                                                         
//This class is used to get some informations of each layers               
//You can get                                                              
//	       * Number of wires,                                           
//             * Interval angle of sense wires,                             
//             * Radius at end plate,                                       
//             * Radius at z=0,                                             
//             * z position at end plate,                                   
//             * Number of holes sense wire skip when it go                 
//               from end plate to another end plate
//
//Define only constructor
//These information is from 140328-COMET-CDC-WireConfig_tilt10_ToIida.xlsx 
//-----------------------------------------------------------//	           

#include <math.h>
#include "LayerInf140328.hh"

#define PI 3.14159265358979
LayerInf::LayerInf(int l){

    double r0[18] = {52.7600,
	             54.3671,
		     55.9738,
		     57.5801,
		     59.1860,
		     60.7917,
		     62.3970,
		     64.0021,
		     65.6069,
		     67.2115,
		     68.8159,
		     70.4201,
		     72.0241,
		     73.6280,
		     75.1709,
		     76.7756,
		     78.3802,
		     79.9846};            //[cm]
    
    double rEP[18] = {53.00,
		      54.60,
		      56.20,
		      57.80,
		      59.40,
		      61.00,
		      62.60,
		      64.20,
		      65.80,
		      67.40,
		      69.00,
		      70.60,
		      72.20,
		      73.80,
		      75.40,
		      77.00,
		      78.60,
		      80.20};
    
    int numberofwire[18] = {396/2,
	                    408/2,
			    420/2,
			    432/2,
			    444/2,
			    456/2,
			    468/2,
			    480/2,
			    492/2,
			    504/2,
			    516/2,
			    528/2,
			    540/2,
			    552/2,
			    564/2,
			    576/2,
			    588/2,
			    600/2};

    double interval[18];                
    for(int i=0;i<18;i++){
	interval[i] = 2*PI/numberofwire[i];
    }

    double zEP[18] = {148.1688/2.,
	              148.7331/2.,
		      149.2973/2.,
		      149.8616/2.,
		      150.4258/2.,
		      150.9900/2.,
		      151.5543/2.,
		      152.1185/2.,
		      152.6828/2.,
		      153.2470/2.,
		      153.8113/2.,
		      154.3755/2.,
		      154.9398/2.,
		      155.5040/2.,
		      156.0683/2.,
		      156.6325/2.,
		      157.1968/2.,
		      157.7610/2.};       //[cm]

//This is used to calculate coordinates of the back of EP

    int NOfskippinghole[18] = {-6,      
	                        6,
			       -6,
			        6,
			       -6,
			        6,
                               -6,
			        6,
                               -6,
			        6,
			       -6,
			        6,
			       -6,
			        6,
                               -7,
			        7,
			       -7,
                                7,};

    NOfWire = numberofwire[l];
    Interval = interval[l];
    R_0 = r0[l];
    R_EP = rEP[l];
    z_EP = zEP[l];
    NOfSkip = NOfskippinghole[l];
    
}
