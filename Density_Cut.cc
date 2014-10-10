//-----------------------------------------------------------//        
//Density_Cut.cc
//
//This function is used to density cut after Neural Net                
//                                                                     
//Argument: vector_of_x-coordinate_before_cut,                         
//	    vector_of_y-coordinate_before_cut,                         
//	    vector_of_x-coordinate_after_cut,                          
//	    vector_of_y-coordinate_after_cut,                          
//	    radius,                                                    
//	    how_many_hits_you_want_in_the_radius                       
//-----------------------------------------------------------//	       

#include <math.h>
#include <stdio.h>
#include <iostream>
#include "Density_CUt.hh"

using namespace std;
void Density_Cut(vector<double> *x, vector<double> *y, vector<double> *clean_x, vector<double> *clean_y, double Radius, int NOfHits){

    for(int i=0;i<x->size();i++){                                           
        int count6 = 0;      
        double R_cut2 = Radius; 
        for(int j=0;j<y->size();j++){
            if(i != j){                                                         
        	double d_ij = sqrt((x->at(j)-x->at(i))*(x->at(j)-x->at(i)) + (y->at(j)-y->at(i))*(y->at(j)-y->at(i))); 
                if(d_ij<R_cut2){                                                  
        	    count6++; 
        	} 
            	if(count6>NOfHits){ 
        	    clean_x->push_back(x->at(i)); 
        	    clean_y->push_back(y->at(i)); 
        	    break; 
        	}    
            } 
        } 
    } 
}    
