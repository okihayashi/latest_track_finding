///Calcurate the closest point of two wires//////////// 

#include <math.h>
#include <iostream>
#include "Distance.hh"

using namespace std;
void Distance(Lineinf* inf){

    //AB m 
    double ABm = ((inf->c1)-(inf->a1))*((inf->b1)-(inf->a1)) + ((inf->c2)-(inf->a2))*((inf->b2)-(inf->a2)) + ((inf->c3)-(inf->a3))*((inf->b3)-(inf->a3));
    
    //AB n 
    double ABn = ((inf->c1)-(inf->a1))*((inf->d1)-(inf->c1)) + ((inf->c2)-(inf->a2))*((inf->d2)-(inf->c2)) + ((inf->c3)-(inf->a3))*((inf->d3)-(inf->c3));

    //m n
    double mn = ((inf->b1)-(inf->a1))*((inf->d1)-(inf->c1)) + ((inf->b2)-(inf->a2))*((inf->d2)-(inf->c2)) + ((inf->b3)-(inf->a3))*((inf->d3)-(inf->c3));

    //m m
    double mm = ((inf->b1)-(inf->a1))*((inf->b1)-(inf->a1)) + ((inf->b2)-(inf->a2))*((inf->b2)-(inf->a2)) + ((inf->b3)-(inf->a3))*((inf->b3)-(inf->a3));

    //n n
    double nn = ((inf->d1)-(inf->c1))*((inf->d1)-(inf->c1)) + ((inf->d2)-(inf->c2))*((inf->d2)-(inf->c2)) + ((inf->d3)-(inf->c3))*((inf->d3)-(inf->c3));

    //Nearest z1
    double z1 = inf->a3 + (ABm * nn - ABn * mn)/(mm * nn - mn * mn) * ((inf->b3) - (inf->a3));
    
    //Nearest z2
    double z2 = inf->c3 + (ABm * mn - ABn * mm)/(nn * mm - mn * mn) * ((inf->d3) - (inf->c3));

    //Return Average of z1 and z2
    inf->z = (z1 + z2)/2.;

    // coordinate of 1st line where z is average of z1 and z2
    inf->x1 = (inf->a1) + ((inf->z) - (inf->a3))/((inf->b3) - (inf->a3)) * ((inf->b1) - (inf->a1));
    inf->y1 = (inf->a2) + ((inf->z) - (inf->a3))/((inf->b3) - (inf->a3)) * ((inf->b2) - (inf->a2));

    // coordinate of 2nd line where z is average of z1 and z2 
    inf->x2 = (inf->c1) + ((inf->z) - (inf->c3))/((inf->d3) - (inf->c3)) * ((inf->d1) - (inf->c1));
    inf->y2 = (inf->c2) + ((inf->z) - (inf->c3))/((inf->d3) - (inf->c3)) * ((inf->d2) - (inf->c2));
    
}    

