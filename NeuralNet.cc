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
#include <TH1D.h>
#include "NeuralNet.hh"

using namespace std;
void NeuralNet(vector<double> *x, vector<double> *y, vector<double> *signalX, vector<double> *signalY, NNParameter* param){

    int k,l;                                                              
    const int NOfHit_cut = x->size();
    const int NOfStep = param->NumberofStep;;                                  
    double theta_tr = 4.90416971470373081e-02*2.;                       
    const double Na = NOfHit_cut;                                               
    const double lambda    = param->lambda1;                                      
    const double kl        = param->a;                                                
    const double ln        = param->b;                                                
    const double alpha     = param->alpha1;                                        
    const double beta      = param->beta1;                                          
    const double C         = param->C1;                                                
    const double T         = param->Temperature;                                       
    const double threshold = param->threshold;                          
    const double R_cut     = param->distance_cut;           
    const double angle_cut = param->A_cut;              
    double T_klnV_ln = 0; double V_kn = 0; double V_ml = 0; double V_mn = 0;                  
    double U_klnV_ln = 0;
    double V[1000][1000] = {{}};
    double R_sum = 0;
    double R_0 = 0;
    int count_segment = 0;
    double U_kln = 0;
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
	if(k != l){
	    double d_kl = sqrt((x->at(l-1)-x->at(k-1))*(x->at(l-1)-x->at(k-1)) + (y->at(l-1)-y->at(k-1))*(y->at(l-1)-y->at(k-1)));  
	    if(d_kl<R_cut){
		//--- select one more hit to calculate T_kln(cost term)
		for(int n=1;n<Na+1;n++){
		    if(n == l || k == n){
			continue;
		    }
		    double d_ln = sqrt((x->at(n-1)-x->at(l-1))*(x->at(n-1)-x->at(l-1)) + (y->at(n-1)-y->at(l-1))*(y->at(n-1)-y->at(l-1))); 
		    if(d_ln<R_cut){
			double theta_kl = acos((x->at(l-1)-x->at(k-1))/d_kl);
			double theta_ln = acos((x->at(n-1)-x->at(l-1))/d_ln);
			double theta_kln = theta_ln - theta_kl;
			if(fabs(theta_kln)<angle_cut){
			    double T_kln = pow(cos(fabs(theta_kln)-theta_tr),lambda)/(pow(d_kl,kl)+pow(d_ln,ln)); 
			    T_klnV_ln += T_kln * V[l][n];
                        }
                    }
		}

		//--- calculate alpha term(constraint term) 
		for(int m=1;m<Na+1;m++){                                                                                                    
		    double d_lm = sqrt((x->at(m-1)-x->at(l-1))*(x->at(m-1)-x->at(l-1)) + (y->at(m-1)-y->at(l-1))*(y->at(m-1)-y->at(l-1)));  
		    if(m != l && d_lm<R_cut){                                                                                               
			V_ml += V[m][l];                                                                                                    
		    }                                                                                                                       
		}                                                                                                                           
		for(int n=1;n<Na+1;n++){
		    double d_nk = sqrt((x->at(n-1)-x->at(k-1))*(x->at(n-1)-x->at(k-1)) + (y->at(n-1)-y->at(k-1))*(y->at(n-1)-y->at(k-1))); 
		    if(n != l && d_nk<R_cut){
			V_kn += V[k][n];
		    }
		}

                ////--- calculate beta term which related to vertex constraint (added in 12th Nov. 2014)
                //for(int n=1;n<Na+1;n++){                                                                                                            
                //    if(n == l || k == n){                                                                                                           
                //	continue;                                                                                                                   
                //    }                                                                                                                               
                //    double d_ln = sqrt((x->at(n-1)-x->at(l-1))*(x->at(n-1)-x->at(l-1)) + (y->at(n-1)-y->at(l-1))*(y->at(n-1)-y->at(l-1)));          
                //    if(d_ln<R_cut){                                                                                                                 
                //	double theta_kl = acos((x->at(l-1)-x->at(k-1))/d_kl);                                                                       
                //	double theta_ln = acos((x->at(n-1)-x->at(l-1))/d_ln);                                                                       
                //	double theta_kln = theta_ln - theta_kl;                                                                                     
                //	if(fabs(theta_kln)<angle_cut){                                                                                              
                //	    //--- Calculate U_kln                                                                                                   
                //	    double A = (x->at(l-1)-x->at(k-1))/(y->at(l-1)-y->at(k-1));                                                           
                //	    double B = (x->at(l-1)-x->at(k-1))*(x->at(k-1)+x->at(l-1))/(2*(y->at(l-1)-y->at(k-1))) + (y->at(k-1)+y->at(l-1))/2.;  
                //	    double C = (x->at(n-1)-x->at(l-1))/(y->at(n-1)-y->at(l-1));                                                           
                //	    double D = (x->at(n-1)-x->at(l-1))*(x->at(l-1)+x->at(n-1))/(2*(y->at(n-1)-y->at(l-1))) + (y->at(l-1)+y->at(n-1))/2.;  
                //	    double x_center = (B-D)/(A-C);                                                                                        
                //	    double y_center = -A*(B-D)/(A-C) + B;                                                                                 
                //	    double R_kln = sqrt((x->at(k-1)-x_center)*(x->at(k-1)-x_center) + (y->at(k-1)-y_center)*(y->at(k-1)-y_center));       
                //	    double R_0 = sqrt(x_center*x_center + y_center*y_center);                                                             
                //	    double R_term = 0;                                                                                                                        
                //	    //if(fabs(R_0-R_kln)<1.){                                                                                              
                //	      R_term = fabs(R_0-R_kln);                                                                                      
                //	    //}                                                                                                                     
                //	    //if(fabs(R_0-R_kln)>1.){                                                                                              
                //	    //    //--- insert some huge number to R_term                                                                           
                //	    //    R_term = DBL_MAX;                                                                                                 
                //	    //                                                                                                                      
                //	    //}                                                                                                                     
                //	                                                                                                                            
                //	    //cout << "R_term = " << R_term << endl;                                                                                
                //	    //cout << "---------------------" << endl;                                                                              
                //	                          
                //	    //double T_kln = pow(cos(fabs(theta_kln)-theta_tr),lambda)/(pow(d_kl,kl)+pow(d_ln,ln));                                   
                //	    U_kln += R_term;
                //        }                                                                                                                           
                //    }                                                                                                                               
                //}                                                                                                                                   

		////--- Calculate new V_kl and reset 
		//if(i<NOfStep-1){
		//    //double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_klnV_ln)); 
		//    double V_kl = tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_kln);  
                //    if(V[k][l] < V_kl){
		//	V[k][l] = V_kl;
		//    }	
		//}
		//if(i == NOfStep-1){
                //    //double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_klnV_ln)); 
		//    double V_kl = tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_kln);
                //    V[k][l] = V_kl;
		//}
		
                if(i < NOfStep*(3/4)){                                                                                         
                    //double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_klnV_ln));  
                    double V_kl = tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml));                      
                    if(V[k][l] < V_kl){                                                                                  
                	V[k][l] = V_kl;                                                                                    
                    }	                                                                                                   
                }                                                                                                        
                
                //--- Calculate R_0
                if(i == NOfStep*(3/4)){
                    cout << "i = " << i << endl;
                    for(int ii=0;ii<Na+1;ii++){
                        for(int jj=0;jj<Na+1;jj++){
                            for(int nn=0;nn<Na+1;nn++){
                                if(ii != jj && jj != nn && nn != ii && V[ii][jj]>0.8 && V[jj][ii]>0.8 && V[jj][nn]>0.8 && V[nn][jj]>0.8){
                                    count_segment++;
                                    double A = (x->at(jj-1)-x->at(ii-1))/(y->at(jj-1)-y->at(ii-1));  
                                    double B = (x->at(jj-1)-x->at(ii-1))*(x->at(ii-1)+x->at(jj-1))/(2*(y->at(jj-1)-y->at(ii-1))) + (y->at(ii-1)+y->at(jj-1))/2.; 
                                    double C = (x->at(nn-1)-x->at(jj-1))/(y->at(nn-1)-y->at(jj-1));  
                                    double D = (x->at(nn-1)-x->at(jj-1))*(x->at(jj-1)+x->at(nn-1))/(2*(y->at(nn-1)-y->at(jj-1))) + (y->at(jj-1)+y->at(nn-1))/2.; 
                                    double x_center = (B-D)/(A-C); 
                                    double y_center = -A*(B-D)/(A-C) + B; 
                                    double R_0 = sqrt(x_center*x_center + y_center*y_center); 

                                    R_sum += R_0;
                                }
                            }    
                        }
                    }
                    R_0 = R_sum/count_segment;
                    cout << "R_0 = " << R_0 << endl;
                }
                if(i > NOfStep*(3/4)){
                    //--- calculate beta term which related to vertex constraint (added in 12th Nov. 2014)                                                                
                    for(int n=1;n<Na+1;n++){                                                                                                                              
                        if(n == l || k == n){                                                                                                                             
                    	continue;                                                                                                                                       
                        }                                                                                                                                                 
                        double d_ln = sqrt((x->at(n-1)-x->at(l-1))*(x->at(n-1)-x->at(l-1)) + (y->at(n-1)-y->at(l-1))*(y->at(n-1)-y->at(l-1)));                            
                        if(d_ln<R_cut){                                                                                                                                   
                    	double theta_kl = acos((x->at(l-1)-x->at(k-1))/d_kl);                                                                                           
                    	double theta_ln = acos((x->at(n-1)-x->at(l-1))/d_ln);                                                                                           
                    	double theta_kln = theta_ln - theta_kl;                                                                                                         
                    	if(fabs(theta_kln)<angle_cut){                                                                                                                  
                    	    //--- Calculate U_kln                                                                                                                       
                    	    double A = (x->at(l-1)-x->at(k-1))/(y->at(l-1)-y->at(k-1));                                                                                 
                    	    double B = (x->at(l-1)-x->at(k-1))*(x->at(k-1)+x->at(l-1))/(2*(y->at(l-1)-y->at(k-1))) + (y->at(k-1)+y->at(l-1))/2.;                        
                    	    double C = (x->at(n-1)-x->at(l-1))/(y->at(n-1)-y->at(l-1));                                                                                 
                    	    double D = (x->at(n-1)-x->at(l-1))*(x->at(l-1)+x->at(n-1))/(2*(y->at(n-1)-y->at(l-1))) + (y->at(l-1)+y->at(n-1))/2.;                        
                    	    double x_center = (B-D)/(A-C);                                                                                                              
                    	    double y_center = -A*(B-D)/(A-C) + B;                                                                                                       
                    	    double R_kln = sqrt((x->at(k-1)-x_center)*(x->at(k-1)-x_center) + (y->at(k-1)-y_center)*(y->at(k-1)-y_center));                             
                    	    double R_term = 0;                                                                                                                          
                    	    //if(fabs(R_0-R_kln)<1.){                                                                                                                   
                    	      R_term = fabs(R_0-R_kln);                                                                                                                 
                    	    //}                                                                                                                                         
                    	    //if(fabs(R_0-R_kln)>1.){                                                                                                                   
                    	    //    //--- insert some huge number to R_term                                                                                               
                    	    //    R_term = DBL_MAX;                                                                                                                    
                    	    //                                                                                                                                         
                    	    //}                                                                                                                                        
                    	                                                                                                                                               
                    	    //cout << "R_term = " << R_term << endl;                                                                                                   
                    	    //cout << "---------------------" << endl;                                                                                                 
                    	                                                                                                                                               
                    	    //double T_kln = pow(cos(fabs(theta_kln)-theta_tr),lambda)/(pow(d_kl,kl)+pow(d_ln,ln));                                                    
                    	    U_kln += R_term;                                                                                                                           
                            }                                                                                                                                            
                        }                                                                                                                                                
                    }                                                                                                                                                    
                    //double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_klnV_ln));         
                    double V_kl = tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_kln);
                    if(V[k][l] < V_kl){
                        V[k][l] = V_kl;                                                                                         
                    }    
                }                                                                                                               
                //if(i == NOfStep-1){                                                                                      
                //    //double V_kl = 1./2. * (1 + tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_klnV_ln));  
                //    double V_kl = tanh((C/T)*T_klnV_ln - (alpha/T)*(V_kn + V_ml) - (beta/T)*U_kln);                      
                //    V[k][l] = V_kl;                                                                                      
                //}                                                                                                        

                T_klnV_ln = 0;
		V_kn = 0;
		V_ml = 0;
		V_mn = 0;
                U_kln = 0;
	    }
	}
    }

    //--- Insert hits that have value over threshold to vector array
    //--- after checking double count 

    int check_i[2000] = {};
    int check_j[2000] = {};

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
        		    if(fabs(x->at(i-1)-signalX->at(k))<DBL_EPSILON){
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
        		    if(fabs(x->at(j-1) == signalX->at(l))<DBL_EPSILON){
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

