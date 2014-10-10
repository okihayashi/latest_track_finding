//-----------------------------------------------------------//              
//Neural Network.hh                                                          
//                                                                           
//This function is used to implement neural network to                       
//track finding algorithm                                                    
//                                                                           
//Argument: vector of x-coordinate before neural net,                        
//	    vector of y-coordinate before neural net,                        
//	    vector of x-coordinate after neural net,                         
//	    ventor of y-coordinate after neural net,                         
//	    structure of parameters
//
//Parameters: NOfStep, lambda, a, b, alpha, beta, C, T, 
//	      V_ij_threshold, distance_cut, angle_cut 
//-----------------------------------------------------------//	             

#ifndef NEURALNET_HH
#define NEURALNET_HH

#include <iostream>
#include <vector>
using namespace std;

struct NNParameter{
    int NumberofStep; double lambda1; double a; double b; double alpha1; double beta1;
    double C1; double Temperature; double threshold; double distance_cut; double A_cut;
};	

void NeuralNet(vector<double> *x, vector<double> *y, vector<double> *signalX, vector<double> *signalY, NNParameter* param);

#endif
