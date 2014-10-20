//-----------------------------------------------------------//           
//Neural Network.cc                                                         
//                                                                        
//This function is used to implement neural network to
//track finding algorithm
//
//Argument: vector of x-coordinate before neural net,
//	    vector of y-coordinate before neural net,
//	    vector of x-coordinate after neural net,
//	    ventor of y-coordinate after neural net,
//	    structure of parameters
//-----------------------------------------------------------//	          

#include <math.h>
#include <TRandom.h>
#include <float.h>
#include "NeuralNet.hh"

using namespace std;
void NeuralNet(vector<double> *x, vector<double> *y, vector<double> *signalX, vector<double> *signalY, NNParameter* param){

    int NOfHit_cut = x->size();

    int k,l;                                                              
    int NOfStep = param->NumberofStep;;                                  
    double theta_tr = 4.90416971470373081e-02*2.;                       
    double Na = NOfHit_cut;                                               
    
    double lambda    = param->lambda1;                                      
    double kl        = param->a;                                                
    double ln        = param->b;                                                
    double alpha     = param->alpha1;                                        
    double beta      = param->beta1;                                          
    double C         = param->C1;                                                
    double T         = param->Temperature;                                       
    double threshold = param->threshold;                          
    double R_cut     = param->distance_cut;           
    double angle_cut = param->A_cut;              

    double T_klnV_ln = 0, V_kn = 0, V_ml = 0, V_mn = 0;                  
    double V[1000][1000] = {{}};

    
    
    //-- Set random value to all V_xy through distance cut -----------//
    
    for(int i=1;i<NOfHit_cut+1;i++){
	for(int j=1;j<NOfHit_cut+1;j++){
	    
	    double d_ij = sqrt((x->at(j-1)-x->at(i-1))*(x->at(j-1)-x->at(i-1)) + (y->at(j-1)-y->at(i-1))*(y->at(j-1)-y->at(i-1)));

	    if(i != j && d_ij<R_cut){
		V[i][j] = gRandom->Uniform(0.,1.);
	    }else{
		V[i][j] = 0;
	    }
	}
    }
    
    //----------------------------------------------------------------//

    //-- Calculate and reset V_xy NOfStep times ----------------------//
    
    for(int i=0;i<NOfStep;i++){
	
	//--- select two hits randomly 
	k = gRandom->Uniform(1., Na-0.6);
	l = gRandom->Uniform(1., Na-0.6); 

	if(k == l && k > l){
	    continue;
	}
	if(k < l){

	    double d_kl = sqrt((x->at(l-1)-x->at(k-1))*(x->at(l-1)-x->at(k-1)) + (y->at(l-1)-y->at(k-1))*(y->at(l-1)-y->at(k-1)));  

	    if(d_kl<R_cut){
		
		//--- select one more hit to calculate T_kln(cost term)
		for(int n=1;n<Na+1;n++){
		    if(n == l || k == n){
			continue;
		    }
                    if(l < n){

			double d_ln = sqrt((x->at(n-1)-x->at(l-1))*(x->at(n-1)-x->at(l-1)) + (y->at(n-1)-y->at(l-1))*(y->at(n-1)-y->at(l-1))); 

			if(d_ln<R_cut){

			    double theta_kl = acos((x->at(l-1)-x->at(k-1))/d_kl);
			    double theta_ln = acos((x->at(n-1)-x->at(l-1))/d_ln);
			    double theta_kln = theta_ln - theta_kl;
			    
			    if(fabs(theta_kln)<angle_cut){
				//--- Calculate T_kln
				double A = (x->at(l-1)-x->at(k-1))/(y->at(l-1)-y->at(k-1));
				double B = (x->at(l-1)-x->at(k-1))*(x->at(k-1)+x->at(l-1))/(2*(y->at(l-1)-y->at(k-1))) + (y->at(k-1)+y->at(l-1))/2.;
				double C = (x->at(n-1)-x->at(l-1))/(y->at(n-1)-y->at(l-1)); 
				double D = (x->at(n-1)-x->at(l-1))*(x->at(l-1)+x->at(n-1))/(2*(y->at(n-1)-y->at(l-1))) + (y->at(l-1)+y->at(n-1))/2.; 
				double x_center = (B-D)/(A-C);
				double y_center = -A*(B-D)/(A-C) + B;
				double R_kln = sqrt((x->at(k-1)-x_center)*(x->at(k-1)-x_center) + (y->at(k-1)-y_center)*(y->at(k-1)-y_center));
				double R_0 = sqrt(x_center*x_center + y_center*y_center);
				double R_term = 0;
				//cout << "xk xl xn = " << x->at(k-1) << " " << x->at(l-1) <<  " " << x->at(n-1) << " " << endl;
				//cout << "yk yl yn = " << y->at(k-1) << " " << y->at(l-1) <<  " " << y->at(n-1) << " " << endl; 
				//cout << "x-center = " << x_center << " , y_center = " << y_center << endl;
				//cout << "R_kln = " << R_kln << endl;
				//cout << "R_0 = " << R_0 << endl; 
				//cout << "R_0-R_kln = " << fabs(R_0-R_kln) << endl;
				
				if(fabs(R_0-R_kln)>0.1){
				    R_term = 1./fabs(R_0-R_kln);
				}
				if(fabs(R_0-R_kln)<0.1){
				    //--- insert some huge number to R_term 
				    R_term = DBL_MAX;
				    
				}
				
				//cout << "R_term = " << R_term << endl;
				//cout << "---------------------" << endl;
				//double T_kln = pow(cos(fabs(theta_kln)-theta_tr),lambda)/(pow(d_kl,kl)+pow(d_ln,ln)) + R_term;
				double T_kln = pow(cos(fabs(theta_kln)-theta_tr),lambda)/(pow(d_kl,kl)+pow(d_ln,ln)); 
				T_klnV_ln += T_kln * V[l][n];

			    }else{
				continue;
			    }

			}else{
			    continue;
			}
		    }
		}

		//--- calculate alpha term(constraint term) 
		for(int m=1;m<Na+1;m++){                                                                                                    
		                                                                                                                            
		    double d_lm = sqrt((x->at(m-1)-x->at(l-1))*(x->at(m-1)-x->at(l-1)) + (y->at(m-1)-y->at(l-1))*(y->at(m-1)-y->at(l-1)));  
		                                                                                                                            
		    if(m != l && d_lm<R_cut){                                                                                               
			V_ml += V[m][l];                                                                                                    
		    }else{                                                                                                                  
			continue;                                                                                                           
		    }                                                                                                                       
		}                                                                                                                           
		
		for(int n=1;n<Na+1;n++){

		    double d_nk = sqrt((x->at(n-1)-x->at(k-1))*(x->at(n-1)-x->at(k-1)) + (y->at(n-1)-y->at(k-1))*(y->at(n-1)-y->at(k-1))); 

		    if(n != l && d_nk<R_cut){
			V_kn += V[k][n];
		    }else{
			continue;
		    }
		}

		//--- Calculate new V_kl and reset 
		if(i<NOfStep-1){
		
		    double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T))); 
		    if(V[k][l] < V_kl){
			V[k][l] = V_kl;
		    }	
		}
		if(i == NOfStep-1){
                    double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T))); 
		    V[k][l] = V_kl;
		}
		
		T_klnV_ln = 0;
		V_kn = 0;
		V_ml = 0;
		V_mn = 0;

	    }else{
		continue;
	    }
	}
    }

    //--- Insert hits that have value over threshold to vector array
    //--- after checking double count 

    int check_i[1000] = {};
    int check_j[1000] = {};

    for(int i=1;i<Na+1;i++){
	for(int j=1;j<Na+1;j++){
	    if(V[i][j] > threshold){

		if(check_i[i] == 0 && check_j[i] == 0){
		    if(!signalX->size()){
                        signalX->push_back(x->at(i-1)); 
			signalY->push_back(y->at(i-1)); 
			check_i[i] = 1;
			check_j[i] = 1;     
		    }			
		    if(signalX->size()){
			int flag1 = 0;
			for(int k=0;k<signalX->size();k++){
			    if(x->at(i-1) == signalX->at(k)){
				flag1 = 1;
			    }
			}
			if(flag1 == 0){
			    signalX->push_back(x->at(i-1));
			    signalY->push_back(y->at(i-1));
			    check_i[i] = 1;
			    check_j[i] = 1;
			}
		    }
		}
		if(check_j[j] == 0 && check_i[j] == 0){    
		    if(!signalX->size()){
                        signalX->push_back(x->at(j-1));   
                        signalY->push_back(y->at(j-1));   
			check_j[j] = 1;
			check_i[j] = 1;
		    }			
		    if(signalX->size()){
			int flag2 = 0;
			for(int l=0;l<signalX->size();l++){
			    if(x->at(j-1) == signalX->at(l)){
				flag2 = 1;
			    }
			}
			if(flag2 == 0){			
			    signalX->push_back(x->at(j-1));  
			    signalY->push_back(y->at(j-1)); 
			    check_j[j] = 1;
			    check_i[j] = 1;
			}
		    }
		}
	    }
	}
    }

}

