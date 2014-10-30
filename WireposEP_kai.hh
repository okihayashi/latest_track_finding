//-----------------------------------------------------------//                 
//WireposEP.hh                                                                  
//                                                                              
//This function is used to get sense wires position at                          
//end plete by using LayerInf140104.cc                                                                    
//                                                                              
//Argument: layer number of sense wire,                                         
//	    cell number of sense wire,                                          
//	    x-coordinate of position of sense wire,                             
//	    y-coordinate of position of sense wire,                             
//	    (z-coordinate of position of sense wire)                            
//-----------------------------------------------------------//	                

#ifndef WIREPOSEP_HH
#define WIREPOSEP_HH

void WireposEP(int layer, double theta, double *x, double *y, double *z);

#endif
