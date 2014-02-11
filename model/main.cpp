//
//  main.cpp
//  Mamcat Test
//
//  Created by Solomon Russell on 12/12/13.
//  Copyright (c) 2013 Solomon Russell. All rights reserved.
//

#include <iostream>
#include "mamcatmodel.h"

/*
int main(int argc, const char * argv[])
{

    // insert code here...
    
    MamcatModel mp;
    
    //mp.
    mp.MType(Mam);
    mp.OutputUnit(Mass);
    mp.PoolNum (3);
    mp.Dose (2.4);
    mp.SSConc1 (3.5);
    mp.BodyWeight(Gram);
    
    mp.MassUnit(M_PDOSE);
    mp.TimeUnit(T_SEC);
    mp.VolUnit(V_ML);
    mp.EndoMU(ENDO_M_G);
*/    
    /* This is a hack, it doesn't work for CAT when
     the pool number is 1, it crashes the program.
     Looks like there's some problem with the mamcat
     code at this point.  Since the answers should be
     the same with 1 compartment just change it to mam
     */
/*    if ((mp.PoolNum() ==1 ) && (mp.MType() == Cat)) {
   		mp.MType(Mam);
    }
    

    
    // loading coefficients and exponentials
    int i;
    for (i = 0; i < mp.PoolNum(); i++) {
        mp.A(i+1, 2);
        mp.L(i+1, -.05 - (i*.001));
    }
    std::cout << mp.Calculate();

    std::cout << "Hello, World!\n";
    return 0;
}
*/

