//....................................................................
// Description of Coulomb blockade
// By using hyperbolic function method from SET-FET model
// For quantum transport modeling in 55nm MOSFETs at 1.2 K
//
// Written by Mr. Yuefeng Chen and Mr. Yuanke Zhang
// Copyright (c) belongs to Mr. Yuanke Zhang and Mr. Yuefeng Chen
//
//  Last edit: 2021.12.23
//....................................................................


#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <thread> 
#include <chrono> 
#include <omp.h>



//#include <time.h>

// #include <stdio.h>
// #include <math.h>
#include "set_final.h"
#include "set_final.cpp"

using namespace std;





    double temp = 1.2;
    double eV = 1.60217733e-19;
    double kB=1.380658e-23;
    int number=13;

    double Ids=0;

    double Cd=8.26738e-19;   
    double Cs=3.9584e-17;
    double Cg=8.36957e-18;
    double Cb=0e-18;
    double Rd=58903.5;
    double Rs=36684.8;
    double C0=0;
    double Q0=0;
    double Vpeak=0.479942;
    double k=100;




int main(){

    double Vd =0.006;
    double Vg = 0.4704;
//    double Vpeak[3]={0,0,0};
    double Ipeak[3]={4.2076E-09, 1.82005e-08, 2.53746e-08};
    double Cd[3]={6.0419e-18, 1e-17, 1e-17};
    double Ctot[3]={1.6961e-17, 1.6961e-17, 1.6961e-17};
    double Beta[3]={77.1477, 77.1477, 77.1477};

            Ids = set_final(Vd,Vg,Ipeak,Cd,Ctot,Beta);
            printf("%0.6e\t\n", Ids);

    return 0;

    
}








