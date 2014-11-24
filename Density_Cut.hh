//-----------------------------------------------------------//
//Density_Cut.hh
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

#ifndef DENSITY_CUT_HH
#define DENSITY_CUT_HH
#include <vector>
#include <iostream>

using namespace std;
void Density_Cut(vector<double> *x, vector<double> *y, vector<double> *clean_x, vector<double> *clean_y, double Radius, int NOfhits);

#endif
