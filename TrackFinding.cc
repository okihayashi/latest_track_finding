#include <iostream>
#include <math.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TH1.h>
#include <TVector2.h>
#include "Wirepos0.hh"
#include "WireposEP_kai.hh"
#include "WireposEP2_kai.hh"
#include "WireposReverse.hh"
#include "NeuralNet.hh"
#include "Density_Cut.hh"
#include "Distance.hh"
#include "LayerInf140328.hh"
#include "DrawDetector.hh"

using namespace std;
int main(int argc, char** argv){

    gROOT->Reset();
    //gROOT->SetBatch();
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gRandom->SetSeed( time(NULL) );

    TApplication app("app", &argc, argv);
    TFile* file = new TFile("/Users/hayashi_oki/Workspace/Track-finding/Helix/NN/Data/signal.140905M02.noise-3.root");
    TTree* t = (TTree*)file->Get("tree");

    //TCanvas* c1 = new TCanvas("c1","c1",10,10,800,800);               
    TCanvas* c2 = new TCanvas("c2","c2",10,10,800,800);
    TCanvas* c3 = new TCanvas("c3","c3",1000,10,1200,600); 
    //TCanvas* c4 = new TCanvas("c4","c4",10,500,800,800);  
    //TCanvas* c5 = new TCanvas("c5","c5",100,10,800,800);  
    //TCanvas* c6 = new TCanvas("c6","c6",10,50,800,800);   
    TCanvas* c7 = new TCanvas("c7","c7",30,40,800,600);
    TCanvas* c8 = new TCanvas("c8","c8",40,40,800,600); 

    //TH1D* h1 = new TH1D("h1","h1",30,-80,80);
    
    int CDCcell_nHits;
    vector<double>* CDCcell_x = 0;
    vector<double>* CDCcell_y = 0;
    vector<double>* CDCcell_z = 0;
    vector<int>* CDCcell_layerID = 0;
    vector<int>* CDCcell_cellID = 0;
    vector<int>* CDCcell_hittype = 0;
    vector<double>* CDCcell_edep = 0;
    vector<double>* CDCcell_px = 0;
    vector<double>* CDCcell_py = 0;
    vector<double>* CDCcell_pz = 0;
    t->SetBranchAddress("CdcCell_nHits",&CDCcell_nHits);
    t->SetBranchAddress("CdcCell_x",&CDCcell_x);
    t->SetBranchAddress("CdcCell_y",&CDCcell_y);
    t->SetBranchAddress("CdcCell_z",&CDCcell_z);
    t->SetBranchAddress("CdcCell_layerID",&CDCcell_layerID);
    t->SetBranchAddress("CdcCell_cellID",&CDCcell_cellID);
    t->SetBranchAddress("CdcCell_hittype",&CDCcell_hittype);
    t->SetBranchAddress("CdcCell_edep",&CDCcell_edep);
    t->SetBranchAddress("CdcCell_px",&CDCcell_px);
    t->SetBranchAddress("CdcCell_py",&CDCcell_py);
    t->SetBranchAddress("CdcCell_pz",&CDCcell_pz);
    int count_NOfEvent = 0;
    int count_efficiency = 0; 
    //double findingefficiency[3445] = {};
    double pt[3445] = {};
    double pz[3445] = {};

    for(int a=0;a<1;a++){
    //for(int a=0;a<t->GetEntries();a++){

	cout << "*****  a = " << a << "  *****" << endl;
	t->GetEntry(1);

	//---check events if the signal exist or not
        int check_event = 0;
	
	for(int b=0;b<CDCcell_nHits;b++){
	    if(CDCcell_hittype->at(b) == 0){
		check_event++;
	    }
	}
	if(check_event<15){
	    cout << "---This event do not have any signals---" << endl;
	}
	if(check_event>=15){
	    vector<double> CDCcell_layerID_cut;
            vector<double> CDCcell_cellID_cut;
	    vector<double> CDCcell_hittype_cut;
	    vector<double> CDCcell_px_signal;
	    vector<double> CDCcell_py_signal;
	    vector<double> CDCcell_pz_signal;
	    vector<double> CDCcell_x_signal;
	    vector<double> CDCcell_y_signal;
	    vector<double> CDCcell_z_signal;
	    
	    for(int i=0;i<CDCcell_nHits;i++){
		if(CDCcell_hittype->at(i) == 0){
		    CDCcell_px_signal.push_back(CDCcell_px->at(i));
                    CDCcell_py_signal.push_back(CDCcell_py->at(i));
		    CDCcell_pz_signal.push_back(CDCcell_pz->at(i));
		    CDCcell_x_signal.push_back(CDCcell_x->at(i));
		    CDCcell_y_signal.push_back(CDCcell_y->at(i)); 
		    CDCcell_z_signal.push_back(CDCcell_z->at(i)); 
		}
	    }	

            //--- Calculate pt
	    pt[a] = sqrt(CDCcell_px_signal.at(1)*CDCcell_px_signal.at(1)*1e6 + CDCcell_py_signal.at(1)*CDCcell_py_signal.at(1)*1e6);
            pz[a] = CDCcell_pz_signal.at(1)*1e3;
            cout << "pt = " << pt[a] << " MeV" << endl;
            cout << "pz = " << pz[a] << " MeV" << endl;
	    
	    //--- ADC cut
	    for(int i=0;i<CDCcell_nHits;i++){
		if(CDCcell_edep->at(i)*1e6 > 7.5){  //[keV]
		    continue;
		}else{
		    CDCcell_layerID_cut.push_back(CDCcell_layerID->at(i));
		    CDCcell_cellID_cut.push_back(CDCcell_cellID->at(i)); 
		    CDCcell_hittype_cut.push_back(CDCcell_hittype->at(i));
		}
	    }
        
            CDCcell_nHits = CDCcell_layerID_cut.size();

            //if(pz[a]<40 || pz[a]>=50){
	    //    cout << "---- skip !!! ----" << endl;
	    //}
	    //if(pz[a]>=40 && pz[a]<50){
		cout << "eventNO:" << count_NOfEvent << endl;
		count_NOfEvent++;
		
                //--- first cut (search crossing hit wires in adjacent layer)---------------------------------------------------//
		
		vector<double> x_0;		   //--* x-coordinate of hit wire                                                     
		vector<double> y_0;                //--* y-coordinate of hit wire                                                     
		vector<double> z_0;                //--* z-coordinate of hit wire(not use now)                                                     
		vector<double> x_02;               //--* x-coordinate of hit wire that path through crossing cut 
		vector<double> y_02;               //--* y-coordinate of hit wire that path through crossing cut 
		vector<double> z_02;               //--* z-coordinate of hit wire that path through crossing cut(not use now) 
		vector<double> x_sig;              //--* x-coordinate of hit wire by signal 
		vector<double> y_sig;              //--* y-coordinate of hit wire by signal 
		vector<double> x_bg;               //--* x-coordinate of hit wire by noise 
		vector<double> y_bg;               //--* x-coordinate of hit wire by noise 
                
                for(int i=0;i<CDCcell_nHits;i++){
                    double x01 = 0; double y01 = 0;
                    int layerID_up_cross = 0;
                    double theta_up_cross = 0;
                    Wirepos0(CDCcell_layerID_cut[i],CDCcell_cellID_cut[i],&x01,&y01);
                    if(CDCcell_hittype_cut[i] == 0){
                        x_sig.push_back(x01);
                        y_sig.push_back(y01);
                    }else{
                        x_bg.push_back(x01);
                        y_bg.push_back(y01);
                    }
                    WireposReverse(x01,y01,&layerID_up_cross,&theta_up_cross);
                    for(int j=0;j<CDCcell_nHits;j++){
                        double x02 = 0; double y02 = 0;
                        int layerID_down_cross = 0;
                        double theta_down_cross = 0;
                        Wirepos0(CDCcell_layerID_cut[j],CDCcell_cellID_cut[j],&x02,&y02);
                        WireposReverse(x02,y02,&layerID_down_cross,&theta_down_cross);
                        if(layerID_up_cross == layerID_down_cross+1 || layerID_up_cross == layerID_down_cross-1){
                            if(fabs(theta_down_cross-theta_up_cross)<=0.18){
                                x_02.push_back(x01);
                                y_02.push_back(y01);
                                break;
                            }
                        }
                    }
                }

		//--- apply neural network and density cut ---------------------------------------------------------------------//
		
		vector<double> signalNN_X;         //--* x-coordinate of hits selected by 1st neural network 
		vector<double> signalNN_Y;         //--* y-coordinate of hits selected by 1st neural network 
		vector<double> signalNN_X_2;       //--* x-coordinate of hits selected by 2nd neural network(not use now)  
		vector<double> signalNN_Y_2;       //--* y-coordinate of hits selected by 2nd neural network(not use now)  
	   
		//---NN Parameter is {NOfStep, lambda, a, b, alpha, beta, C, T, V_ij_threshold, distance_cut, angle_cut}
		NNParameter param1 = {500000, 1., 1., 1., 10., 10., 10, 1., 0.8, 6., 0.1};
		NeuralNet(&x_02, &y_02, &signalNN_X, &signalNN_Y, &param1);
		
		vector<double> signalNN_X_cut;     //--* x-coordinate of hits path through the density cut
		vector<double> signalNN_Y_cut;     //--* y-coordinate of hits path through the density cut 
		
		//--- 1st density cut
		if(signalNN_X.size() > 90){
		    Density_Cut(&signalNN_X, &signalNN_Y, &signalNN_X_cut, &signalNN_Y_cut, 20., 15);
		}else{
		    //--- if number of hits is not so large, density cut will be not applied
		    cout << "No more cuts apply" << endl; 
		    signalNN_X_cut = signalNN_X;
		    signalNN_Y_cut = signalNN_Y;
		}

		//--- 2nd neural net(not use now)
		//NNParameter param2 = {50000, 15., 1.0, 1.0, 15., 0.1, 10., 1.0, 0.8, 50., 0.1};
		//NeuralNet(&signalNN_X_cut, &signalNN_Y_cut, &signalNN_X_2, &signalNN_Y_2, &param2);

		vector<double> signalNN_X_2_cut;   //--* x-coordinate of hits path through the 2nd density cut  
		vector<double> signalNN_Y_2_cut;   //--* y-coordinate of hits path through the 2nd density cut  
		
		//--- 2nd density cut
		if(signalNN_X_cut.size() > 70){
		    Density_Cut(&signalNN_X_cut,&signalNN_Y_cut, &signalNN_X_2_cut, &signalNN_Y_2_cut, 40., 30);
		}else{
		    //--- if number of hits is not so large, density cut will be not applied 
		    cout << "No cut is applied in 2nd density cut" << endl;
		    signalNN_X_2_cut = signalNN_X_cut;
		    signalNN_Y_2_cut = signalNN_Y_cut;
		}
                
		//--------------------------------------------------------------------------------------------------------------// 
		
		//--- Introduce 3D information ---------------------------------------------------------------------------------// 
                
                vector<double>   signalNN_X_2_cutEP;
                vector<double>   signalNN_Y_2_cutEP;
	       	vector<double>   signalNN_Z_2_cutEP;
		vector<double>   signalNN_X_2_cutEP2; 
		vector<double>   signalNN_Y_2_cutEP2; 
		vector<double>   signalNN_Z_2_cutEP2; 
		vector<int>      signalNN_layerID;
                vector<double>   signalNN_theta;
		vector<double>   signal_distance_X;
		vector<double>   signal_distance_Y;
		vector<double>   signal_distance_Z;
                vector<double>   signal_distance_X_cut; 
                vector<double>   signal_distance_Y_cut; 
                vector<double>   signal_distance_Z_cut;
                vector<double>   signal_distance_X_cluster;
                vector<double>   signal_distance_Y_cluster;
                vector<double>   signal_distance_Z_cluster; 
                vector<TVector2> signal_distance_origin_X;
                vector<TVector2> signal_distance_origin_Y;
		
                //for(int i=0;i<signalNN_X_2_cut.size();i++){                                                                                                                               
                //    int layerID_re = 0;                                                                                                                                                   
                //    double theta_re = 0;                                                                                                                                                  
                //    double x = 0; double y = 0; double z = 0;                                                                                                                             
                //    double x2 = 0; double y2 = 0; double z2 = 0;                                                                                                                          
                //    WireposReverse(x_sig[i],y_sig[i],&layerID_re,&theta_re);                                                                                        
                //    WireposEP(layerID_re,theta_re,&x,&y,&z);                                                                                                                              
                //    WireposEP2(layerID_re,theta_re,&x2,&y2,&z2);                                                                                                                          
                //    signalNN_layerID.push_back(layerID_re);                                                                                                                               
                //    signalNN_theta.push_back(theta_re);                                                                                                                                   
                //    signalNN_X_2_cutEP.push_back(x);                                                                                                                                      
                //    signalNN_Y_2_cutEP.push_back(y);                                                                                                                                      
                //    signalNN_Z_2_cutEP.push_back(z);                                                                                                                                      
                //    signalNN_X_2_cutEP2.push_back(x2);                                                                                                                                    
                //    signalNN_Y_2_cutEP2.push_back(y2);                                                                                                                                    
                //    signalNN_Z_2_cutEP2.push_back(z2);                                                                                                                                    
                //}                                                                                                                                                                         
                
                //--- transform. (x,y) to (layerID, theta)
		for(int i=0;i<signalNN_X_2_cut.size();i++){
		    int layerID_re = 0;
		    double theta_re = 0;
		    double x = 0; double y = 0; double z = 0;
		    double x2 = 0; double y2 = 0; double z2 = 0;
		    WireposReverse(signalNN_X_2_cut[i],signalNN_Y_2_cut[i],&layerID_re,&theta_re);
		    WireposEP(layerID_re,theta_re,&x,&y,&z);
		    WireposEP2(layerID_re,theta_re,&x2,&y2,&z2);
		    signalNN_layerID.push_back(layerID_re);
		    signalNN_theta.push_back(theta_re);
		    signalNN_X_2_cutEP.push_back(x);
                    signalNN_Y_2_cutEP.push_back(y);
                    signalNN_Z_2_cutEP.push_back(z);
                    signalNN_X_2_cutEP2.push_back(x2); 
                    signalNN_Y_2_cutEP2.push_back(y2); 
                    signalNN_Z_2_cutEP2.push_back(z2);
		}

                int count_ID = 0;
                for(int i=0;i<signalNN_X_2_cut.size();i++){
		    for(int j=0;j<signalNN_X_2_cut.size();j++){
			if((signalNN_layerID[i]+1 == signalNN_layerID[j]) || (signalNN_layerID[i]-1 == signalNN_layerID[j])){
			    if(fabs(signalNN_theta[j]-signalNN_theta[i])<=0.18){    
				double x_near_i = 0; double y_near_i = 0; double x_near_j = 0; 
				double y_near_j = 0; double z_near = 0;
				LayerInf layerinf3(signalNN_layerID[i]);
			        Lineinf inf1 = {signalNN_X_2_cutEP2[i],signalNN_Y_2_cutEP2[i],signalNN_Z_2_cutEP2[i],
						signalNN_X_2_cutEP[i], signalNN_Y_2_cutEP[i], signalNN_Z_2_cutEP[i],
						signalNN_X_2_cutEP2[j],signalNN_Y_2_cutEP2[j],signalNN_Z_2_cutEP2[j],
						signalNN_X_2_cutEP[j], signalNN_Y_2_cutEP[j], signalNN_Z_2_cutEP[j],
						x_near_i, x_near_j, y_near_i, y_near_j, z_near};
				Distance(&inf1);
                                signal_distance_origin_X.push_back(TVector2(signalNN_X_2_cut[i],signalNN_X_2_cut[j]));
                                signal_distance_origin_Y.push_back(TVector2(signalNN_Y_2_cut[i],signalNN_Y_2_cut[j]));
                                signal_distance_X.push_back((x_near_i+x_near_j)/2.);
				signal_distance_Y.push_back((y_near_i+y_near_j)/2.);
				signal_distance_Z.push_back(z_near);
                                cout << "signal_distance_origin_X[" << count_ID << "] = " << signal_distance_origin_X[count_ID].X() << " , " << signal_distance_origin_Y[count_ID].Y() << endl;;
                                count_ID++;
			    }
			}
		    }
		}
                
                //--- Remove hit points using z position
                for(int i=0;i<signal_distance_X.size();i++){
                    int count_distance = 0;
                    for(int j=0;j<signal_distance_X.size();j++){
                        if(fabs(signal_distance_X[i]-signal_distance_X[j])<5. && fabs(signal_distance_Y[i]-signal_distance_Y[j])<5. && fabs(signal_distance_Z[i]-signal_distance_Z[j])<15){
                            count_distance++;
                        }
                        if(count_distance>10){
                            signal_distance_X_cluster.push_back(signal_distance_X[i]);
                            signal_distance_Y_cluster.push_back(signal_distance_Y[i]);
                            signal_distance_Z_cluster.push_back(signal_distance_Z[i]);
                            bool flag_i = 0;
                            bool flag_j = 0;
                            for(int k=0;k<signal_distance_X_cut.size();k++){
                                if(fabs(signal_distance_X_cut[k]-signal_distance_origin_X[i].X())<DBL_EPSILON && fabs(signal_distance_Y_cut[k]-signal_distance_origin_Y[i].X())<DBL_EPSILON){
                                    flag_i = 1;
                                }
                                if(fabs(signal_distance_X_cut[k]-signal_distance_origin_X[i].Y())<DBL_EPSILON && fabs(signal_distance_Y_cut[k]-signal_distance_origin_Y[i].Y())<DBL_EPSILON){
                                    flag_j = 1;
                                }
                            }
                            if(!flag_i){
                                signal_distance_X_cut.push_back(signal_distance_origin_X[i].X());
                                signal_distance_Y_cut.push_back(signal_distance_origin_Y[i].X());
                            }
                            if(!flag_j){
                                signal_distance_X_cut.push_back(signal_distance_origin_X[i].Y());
                                signal_distance_Y_cut.push_back(signal_distance_origin_Y[i].Y());
                            }
                            break;
                        }
                    } 
                }
		

                //--------------------------------------------------------------------------------------------------------------//
		
		//--- count signal and noise hits to calculate efficiency ------------------------------------------------------//

		double countsig = 0.;
		double countsig_2 = 0.;
		double countsig_3 = 0.;
                double countsig_4 = 0.;
                double countsig_5 = 0.;

		cout << "signalsize is " << x_sig.size() << endl;
		
		for(int i=0;i<x_sig.size();i++){                                                                   
                    for(int j=0;j<x_02.size();j++){                                                          
                	if(fabs(x_sig[i]-x_02[j])<DBL_EPSILON && fabs(y_sig[i]-y_02[j])<DBL_EPSILON){  
                	    countsig_5++;                                                                          
                	}                                                                                          
                    }                                                                                              
                }                                                                                                  
                for(int i=0;i<x_sig.size();i++){                                            
		    for(int j=0;j<signalNN_X.size();j++){                               
			if(fabs(x_sig[i]-signalNN_X[j])<DBL_EPSILON && fabs(y_sig[i]-signalNN_Y[j])<DBL_EPSILON){ 
			    countsig_3++;                                                       
			}                                                                   
		    }                                                                       
		}                                                                          
		for(int i=0;i<x_sig.size();i++){
		    for(int j=0;j<signalNN_X_cut.size();j++){
                        if(fabs(x_sig[i]-signalNN_X_cut[j])<DBL_EPSILON && fabs(y_sig[i]-signalNN_Y_cut[j])<DBL_EPSILON){
			    countsig++;
			}
		    }
		}
		for(int i=0;i<x_sig.size();i++){                                             
		    for(int j=0;j<signalNN_X_2_cut.size();j++){                                    
			if(fabs(x_sig[i]-signalNN_X_2_cut[j])<DBL_EPSILON && fabs(y_sig[i]-signalNN_Y_2_cut[j])<DBL_EPSILON){          
			    countsig_2++;                                                        
			}                                                                    
		    }                                                                        
		}                                                                            
                for(int i=0;i<x_sig.size();i++){                                                                                           
                    for(int j=0;j<signalNN_X_2_cut.size();j++){                                                                            
                	if(fabs(x_sig[i]-signal_distance_X_cut[j])<DBL_EPSILON && fabs(y_sig[i]-signal_distance_Y_cut[j])<DBL_EPSILON){           
                	    countsig_4++;                                                                                                    
                	}                                                                                                                  
                    }                                                                                                                      
                }

		double countbg = signalNN_X_cut.size() - countsig;
		double countbg_2 = signalNN_X_2_cut.size() - countsig_2;
		double countbg_3 = signalNN_X.size() - countsig_3;
		double countbg_4 = signal_distance_X_cut.size() - countsig_4;
                double countbg_5 = x_02.size() - countsig_5;

		double efficiency = (countsig/x_sig.size())*100.; 
		double efficiency_2 = (countsig_2/x_sig.size())*100.;
		double efficiency_3 = (countsig_3/x_sig.size())*100.;
		double efficiency_4 = (countsig_4/x_sig.size())*100.;
                double efficiency_5 = (countsig_5/x_sig.size())*100.;
                double noise_reduction = 100. - (countbg/x_bg.size())*100.;
		double noise_reduction_2 = 100. - (countbg_2/x_bg.size())*100.;
		double noise_reduction_3 = 100. - (countbg_3/x_bg.size())*100.;
		double noise_reduction_4 = 100. - (countbg_4/x_bg.size())*100.;
		double noise_reduction_5 = 100. - (countbg_5/x_bg.size())*100.;

		cout << "-------" << endl;
                cout << "After cross cut, Efficiency is" << efficiency_5 << endl;
		cout << "NOfHit: " << x_02.size() << endl;
		cout << "Noise reductoin rate is " << noise_reduction_5 << endl;
                cout << "signal vs noise = " << countsig_5 << ":" << countbg_5 << endl;
		cout << "-------" << endl;
		cout << "-------" << endl;
		cout << "After NN, Efficiency is:" << efficiency_3 << " [%] " << endl;
		cout << "NOfhit is " << signalNN_X.size() << endl;
		cout << "Noise reduction rate is " << noise_reduction_3 << " [%] " << endl;
		cout << "signal vs noise = " << countsig_3 << ":" << countbg_3 << endl; 
                cout << "-------" << endl;
		cout << "--------" << endl; 
		cout << "In 1st density cut, Efficiency is:" << efficiency << " [%] " << endl;
	        cout << "NOfHit is " << signalNN_X_cut.size() << endl;   
		cout << "Noise reduction rate is " << noise_reduction << " [%] " << endl;
	        cout << "signal vs noise = " << countsig << ":" << countbg << endl; 	
                cout << "--------" << endl;
		cout << "--------" << endl;                                             
		cout << "In 2nd density cut, Efficiency is:" << efficiency_2 << " [%] " << endl;  
		cout << "NOfHit is " << signalNN_X_2_cut.size() << endl;    
		cout << "Noise reduction rate is " << noise_reduction_2 << " [%] " << endl;                                  
		cout << "signal vs noise = " << countsig_2 << ":" << countbg_2 << endl; 
                cout << "--------" << endl;
                cout << "--------" << endl;                                                                                   
                cout << "In z cut, Efficiency is:" << efficiency_4 << " [%] " << endl;                                      
                cout << "NOfHit is " << signal_distance_X_cut.size() << endl;                                                      
                cout << "Noise reduction rate is " << noise_reduction_4 << " [%] " << endl;                                   
                cout << "signal vs noise = " << countsig_4 << ":" << countbg_4 << endl; 
                cout << "--------" << endl;                                                                                   


		//--- conditions of 'success'
		if(efficiency_2>=50. && (countsig_2/countbg_2)>=2.){
		    count_efficiency++;
		}
                //if(efficiency_2>=50.){
	        //    count_efficiency++;
		//}
		
		cout << endl;
		cout << "Finally," << endl;
		cout << "number of events 80 [%] over :" << count_efficiency << endl;
		cout << "signal vs noise = " << countsig_4 << ":" << countbg_4 << endl;
		cout << "======================================" << endl;
		
		//findingefficiency[a] = efficiency_2;

		//h1->Fill(pt[a]);
		
		//--------------------------------------------------------------------------------------------------------------// 
		
		//--- clear vectors
		
		//x_0.clear(); 
		//y_0.clear(); 
		//z_0.clear(); 
		//x_02.clear();
		//y_02.clear();
		//z_02.clear();
		//x_sig.clear();
		//y_sig.clear();
		//x_bg.clear(); 
		//y_bg.clear(); 
		//signalNN_X.clear(); 
		//signalNN_Y.clear(); 
		//signalNN_X_2.clear();
		//signalNN_Y_2.clear();
		//signalNN_X_cut.clear();
		//signalNN_Y_cut.clear();
		//signalNN_X_2_cut.clear();
		//signalNN_Y_2_cut.clear();
                //signalNN_X_2_cutEP.clear(); 
                //signalNN_Y_2_cutEP.clear(); 
                //signalNN_Z_2_cutEP.clear(); 
                //signalNN_layerID.clear();      
                //signalNN_theta.clear();     
                //signal_distance_X.clear();  
                //signal_distance_Y.clear();  
                //signal_distance_Z.clear();  
                //signal_distance_X_cut.clear(); 
                //signal_distance_Y_cut.claer(); 
                //signal_distance_Z_cut.clear();
                //signal_distance_X_cluster.clear(); 
                //signal_distance_Y_cluster.clear(); 
                //signal_distance_Z_cluster.clear(); 
                //signalNN_X_2_cutEP2.clear();  
                //signalNN_Y_2_cutEP2.clear();  
                //signalNN_Z_2_cutEP2.clear();  


		//--- drawing --------------------------------------------------------------------------------------------------// 
		
		c2->cd();                                                                        
                c2->Divide(2,2);
                c2->cd(1);
		//TGraph* gsigNN = new TGraph(signalNN_X.size(), &signalNN_X[0], &signalNN_Y[0]);  
		TGraph* gsigNN = new TGraph(x_02.size(), &x_02[0], &y_02[0]);   
		gsigNN->Draw("ap");                                                          
		gsigNN->SetMarkerStyle(7);                                                       
		gsigNN->SetMarkerSize(0.3);                                                      
		gsigNN->SetMarkerColor(4);                                                       
		//gsigNN->SetTitle("after 1st NN");
		gsigNN->SetTitle(0);
		gsigNN->GetXaxis()->SetTitle("X [cm]");
		gsigNN->GetYaxis()->SetTitle("Y [cm]");
		gsigNN->SetMaximum(90);
		gsigNN->SetMinimum(-90);
		gsigNN->GetXaxis()->SetLimits(-90,90);
		DrawDetector();

		c2->cd(2);                                                       
		//TGraph* gcut = new TGraph(x_02.size(), &x_02[0], &y_02[0]);  
		//gcut->Draw("p,same");                                           
		//gcut->SetMarkerStyle(4);                                        
		//gcut->SetMarkerSize(0.5);                                       
		//gcut->SetMarkerColor(4);                                        
	
		//TGraph* gsig = new TGraph(CDCcell_x_signal.size(), &CDCcell_x_signal[0], &CDCcell_y_signal[0]);  
		TGraph* gsig = new TGraph(x_sig.size(), &x_sig[0], &y_sig[0]);
		gsig->Draw("ap");                                           
		gsig->SetMarkerStyle(7);                                        
		gsig->SetMarkerSize(0.3);                                       
		gsig->SetMarkerColor(2);                                        
		//gsig->SetTitle("noise 9%");
		gsig->SetTitle(0);
		gsig->GetXaxis()->SetTitle("X [cm]"); 
		gsig->GetYaxis()->SetTitle("Y [cm]"); 
		gsig->SetMaximum(90);                
		gsig->SetMinimum(-90);               
		gsig->GetXaxis()->SetLimits(-90,90); 

		TGraph* gbg = new TGraph(x_bg.size(), &x_bg[0], &y_bg[0]);      
		gbg->Draw("p,same");                                            
		gbg->SetMarkerStyle(7);                                         
		gbg->SetMarkerSize(0.3);                                        
		gbg->SetMarkerColor(4);                                         
                DrawDetector();

		c2->cd(3);
		TGraph* gsigNN2 = new TGraph(signalNN_X.size(), &signalNN_X[0], &signalNN_Y[0]);   
		gsigNN2->Draw("ap");                                                                 
		gsigNN2->SetMarkerStyle(7);                                                              
		gsigNN2->SetMarkerSize(0.3);                                                             
		gsigNN2->SetMarkerColor(4);                                                              
		gsigNN2->SetTitle("after NN");
		gsigNN2->GetXaxis()->SetTitle("X [cm]");  
		gsigNN2->GetYaxis()->SetTitle("Y [cm]");  
		gsigNN2->SetMaximum(90);                 
		gsigNN2->SetMinimum(-90);                
		gsigNN2->GetXaxis()->SetLimits(-90,90);  
                DrawDetector();

		//c5->cd();                                                                                 
		//TGraph* gsigNN_cut = new TGraph(signalNN_X_cut.size(), &signalNN_X_cut[0], &signalNN_Y_cut[0]);    
		//gsigNN_cut->Draw("ap");                                                                  
		//gsigNN_cut->SetMarkerStyle(4);                                                               
		//gsigNN_cut->SetMarkerSize(0.3);                                                              
		//gsigNN_cut->SetMarkerColor(4);                                                               
		////gsigNN_cut->SetTitle("after 1st NN cut");
		//gsigNN_cut->SetTitle(0);
		//gsigNN_cut->GetXaxis()->SetTitle("X [cm]");  
		//gsigNN_cut->GetYaxis()->SetTitle("Y [cm]");  
		//gsigNN_cut->SetMaximum(90);                  
		//gsigNN_cut->SetMinimum(-90);                 
		//gsigNN_cut->GetXaxis()->SetLimits(-90,90);   
		
		c2->cd(4);                                                                                 
		TGraph* gsigNN2_cut = new TGraph(signalNN_X_2_cut.size(), &signalNN_X_2_cut[0], &signalNN_Y_2_cut[0]);    
		gsigNN2_cut->Draw("ap");                                                                  
		gsigNN2_cut->SetMarkerStyle(7);                                                               
		gsigNN2_cut->SetMarkerSize(0.4);                                                              
		gsigNN2_cut->SetMarkerColor(4);                                                               
		//gsigNN2_cut->SetTitle("after 2nd NN cut: Result");
		gsigNN2_cut->SetTitle(0);
		gsigNN2_cut->GetXaxis()->SetTitle("X [cm]");  
		gsigNN2_cut->GetYaxis()->SetTitle("Y [cm]");  
		gsigNN2_cut->SetMaximum(90);                  
		gsigNN2_cut->SetMinimum(-90);                 
		gsigNN2_cut->GetXaxis()->SetLimits(-90,90);   

                TGraph* gdis = new TGraph(signal_distance_X.size(), &signal_distance_X[0], &signal_distance_Y[0]);
		gdis->Draw("p,same");
		gdis->SetMarkerStyle(3);
		gdis->SetMarkerSize(0.3);
		gdis->SetMarkerColor(6);
                DrawDetector();

		c7->cd();
		TGraph* gyz = new TGraph(signal_distance_X.size(), &signal_distance_Z[0], &signal_distance_Y[0]);
		gyz->Draw("ap");
		gyz->SetMarkerStyle(3);
		gyz->SetMarkerSize(0.4);
		gyz->SetMarkerColor(4);
		gyz->SetMaximum(90);
		gyz->SetMinimum(-90);
		gyz->GetXaxis()->SetLimits(-80,80);
		gyz->GetXaxis()->SetTitle("Z [cm]");
		gyz->GetYaxis()->SetTitle("Y [cm]");

		TGraph* gyz_sig = new TGraph(CDCcell_x_signal.size(), &CDCcell_z_signal[0], &CDCcell_y_signal[0]);
                gyz_sig->Draw("p,same");
		gyz_sig->SetMarkerStyle(7);
		gyz_sig->SetMarkerSize(0.4);
		gyz_sig->SetMarkerColor(2);

                c8->cd();
		TGraph* gsigNN2_Zcluster = new TGraph(signal_distance_X_cluster.size(),&signal_distance_Z_cluster[0],&signal_distance_Y_cluster[0]);
		gsigNN2_Zcluster->Draw("ap");
		gsigNN2_Zcluster->SetMarkerStyle(7);
		gsigNN2_Zcluster->SetMarkerSize(0.4);
		gsigNN2_Zcluster->SetMarkerColor(4);
		gsigNN2_Zcluster->SetMaximum(90);                 
                gsigNN2_Zcluster->SetMinimum(-90);                
                gsigNN2_Zcluster->GetXaxis()->SetLimits(-80,80);  
                gsigNN2_Zcluster->GetXaxis()->SetTitle("Z [cm]"); 
                gsigNN2_Zcluster->GetYaxis()->SetTitle("Y [cm]"); 

                c3->cd();
                c3->Divide(2,1);
                c3->cd(2);                
                TGraph* g_zcut = new TGraph(signal_distance_X_cut.size(), &signal_distance_X_cut[0], &signal_distance_Y_cut[0]);          
                g_zcut->Draw("ap");                                                                      
                g_zcut->SetMarkerStyle(4);                                                               
                g_zcut->SetMarkerSize(0.3);                                                              
                g_zcut->SetMarkerColor(4);                                                               
                g_zcut->SetTitle("after z cut");                                                            
                g_zcut->GetXaxis()->SetTitle("X [cm]");                                                  
                g_zcut->GetYaxis()->SetTitle("Y [cm]");                                                  
                g_zcut->SetMaximum(90);                                                                  
                g_zcut->SetMinimum(-90);                                                                 
                g_zcut->GetXaxis()->SetLimits(-90,90);                                                   
                DrawDetector();     

                c3->cd(1);                                                                                                 
                TGraph* gsigNN2_cut2 = new TGraph(signalNN_X_2_cut.size(), &signalNN_X_2_cut[0], &signalNN_Y_2_cut[0]);     
                gsigNN2_cut2->Draw("ap");                                                                                   
                gsigNN2_cut2->SetMarkerStyle(4);                                                                            
                gsigNN2_cut2->SetMarkerSize(0.4);                                                                           
                gsigNN2_cut2->SetMarkerColor(4);                                                                            
                //gsigNN2_cut2->SetTitle("after 2nd NN cut: Result");                                                       
                gsigNN2_cut2->SetTitle("before z cut");                                                                                   
                gsigNN2_cut2->GetXaxis()->SetTitle("X [cm]");                                                               
                gsigNN2_cut2->GetYaxis()->SetTitle("Y [cm]");                                                               
                gsigNN2_cut2->SetMaximum(90);                                                                               
                gsigNN2_cut2->SetMinimum(-90);                                                                              
                gsigNN2_cut2->GetXaxis()->SetLimits(-90,90);                                                                
                DrawDetector();

    	    //}
	    CDCcell_x_signal.clear();
	    CDCcell_y_signal.clear();
	    CDCcell_z_signal.clear();
	    CDCcell_px_signal.clear();
	    CDCcell_py_signal.clear();
	    CDCcell_pz_signal.clear();
	    CDCcell_layerID_cut.clear(); 
	    CDCcell_cellID_cut.clear();
	    CDCcell_hittype_cut.clear();
	}
    }

    double Result = count_efficiency/double(count_NOfEvent+1);
    double error = sqrt(count_efficiency)/double(count_NOfEvent+1);
    cout << "**----------------------------**" << endl;
    cout << "Result: " << Result << "± " << error << endl;
    cout << "**----------------------------**" << endl; 
    //TGraph* geff = new TGraph(3445, pt, findingefficiency);
    //c1->cd();
    //geff->Draw("ap");
    //geff->SetMarkerStyle(2);
    //geff->SetMarkerSize(1);
    //geff->SetMarkerColor(2);
    //geff->GetXaxis()->SetLimits(0,110);
    //geff->SetMaximum(100);
    //geff->SetMinimum(0);    

    //c7->cd();
    //h1->Draw();
    ////h1->SetFillStyle(1001);
    ////h1->SetFillColor(6);
    ////h1->SetLineColor(6);
    ////h1->SetStats(0);
    //h1->SetTitle(0);
    //h1->GetXaxis()->SetTitle("z [cm]");

    //c1->Update();
    c2->Update();
    c3->Update();
    //c4->Update();
    //c5->Update();
    //c6->Update();
    c7->Update();
    c8->Update();
    app.Run(); } 


    
