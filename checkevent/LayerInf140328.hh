//-----------------------------------------------------------//              
//LayerInf.hh                                                                
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
//These information is from 140104-COMET-CDC-WireConfig_d2_v2.xlsx           
//-----------------------------------------------------------//	             


#ifndef LayerInf_HH
#define LayerInf_HH

class LayerInf {
    private:
	int NOfWire;
	double Interval;
	double R_EP;
	double R_0;
	double z_EP;
        double NOfSkip;
    
    public:
	LayerInf(int l);

	double GetNOfWire()               { return NOfWire;  }
        double GetInterval()              { return Interval; }
	double GetR_EP()                  { return R_EP;     }
        double GetR_0()                   { return R_0;      }
	double Getz_EP()                  { return z_EP;     }
        double GetNOfSkip()               { return NOfSkip;  }
};

#endif
