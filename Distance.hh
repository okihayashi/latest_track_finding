///Calcurate the closest point of two wires////////////

#ifndef __INCLUDE_DISTANCE_HH__
#define __INCLUDE_DISTANCE_HH__
/*-------------------------------------------------------
 * a1,a2,a3 are (x,y,z) of upper wire in the surface
 * b1,b2,b3 are (x,y,z) of upper wire in the back side
 * c1,c2,c3 are         of lower             surface
 * d1,d2,d3 are         of lower             back side
--------------------------------------------------------*/

struct Lineinf{
   double a1; double a2; double a3;  
   double b1; double b2; double b3;  
   double c1; double c2; double c3;  
   double d1; double d2; double d3;  
   double &x1;double &x2;double &y1; 
   double &y2;double &z;            
};    
    
void Distance(Lineinf* inf);

#endif 
