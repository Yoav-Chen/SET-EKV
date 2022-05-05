//..........................................................................
// Cryogenic model based on EKV2.6
// DIBL, DITS related codes are referenced from the EKV3 and BSIM4 model
// Cryogenic mobility correction
//
// Written by Mr. Yuanke Zhang
// Copyright (c) belongs to Mr. Yuanke Zhang
//
//  Last edit: 2021.9.21
//..........................................................................


#ifndef __FUNC_H_
#define __FUNC_H_

#include <iostream>
#include <cstdio>
#include <math.h>
#include <vector>
//#include "set_model.h"
#include "set_model.h"
#include "set_final.h"
#define eps 1e-40
using namespace std;

class BaseFunc
{
    public:
        virtual double func(vector<double>& params, vector<double>& input, int inputSize) = 0;
};

class MatlabFunc : public BaseFunc
{
    public:
        virtual double func(vector<double>& params, vector<double>& input, int inputSize)
        {
            if(input.size() != inputSize){
                printf("input size must be %d.\n", inputSize);
                exit(1);
            }
            if(params.size() != 17){
                printf("parameters size must be right.\n");
                exit(1);
            }
            // input
            double VG = input[0];
            double VD = input[1];
        //    double W = input[2];
        //    double L = input[3];
        
            double VS = 0;
            double VB = 0;
            double Idot=0;
            double temp = 1.2;
            double eV = 1.60217733e-19;
            double kB = 1.380658e-23;
            int number = 13;
            int i;
            double Ipeak[3]={params[5],params[6],params[7]};
            double Cd[3]={params[8],params[9],params[10]};
            double Ctot[3]={params[11],params[12],params[13]};
            double Beta[3] = {params[14],params[15],params[16]};
            

        //     double Cgdot = 8.54836e-18;
        //     double Cddot = 2.53503e-17;
        //     double Csdot = 7.22236e-18;
        // //    double Cbdot = params[8]; 
        //     double Cbdot = 0;
        //     double Rd = 45741.9;
        //     double Rs = 22713.1;


            double C0=0;
            double Q0=0;
            double k=100;
/*             double amp[5]={};
            amp[1]=params[11];
            amp[2]=params[12];
            amp[3]=params[13];
            amp[4]=params[14];
            amp[5]=params[15]; */
            double distance = 0;
            double ref_vg =0;

//            double Vpeak = params[10];
            

            int    i_state[number];
            double P[number];
            double P_sum=0;
            double Vi[number];
            double d_Fdi[number];
            double d_Fid[number];
            double d_Fsi[number];
            double d_Fis[number];

            double Gamma_DR[number];
            double Gamma_DL[number];
            double Gamma_SL[number];
            double Gamma_SR[number];

            



            double Vt = 0.00010341; // double Vt = 0.006600743253@77K; 0.0003626@4.2K;
            double L = 0.055E-6;//double L = 10E-6;
            double W = 0.3E-6; // double W = 10E-6;
            double VTO = params[0]; // double VTO = 0.15382458;
            double GAMMA = params[1]; // double GAMMA = 0.00058673;
            double PHI = params[2]; // double PHI = 1.25796;
            double KP = params[3]; // double KP = 20E-6;
            double THETA = params[4]; // double THETA = 50.0E-3;

            double VGprime = VG - VTO + PHI + GAMMA * sqrt(PHI);
            double VP = VGprime - PHI - GAMMA * (sqrt(VGprime+(GAMMA/2.0)*(GAMMA/2.0))-(GAMMA/2.0));

            // Slope factor
            double n = 1.0 + GAMMA/(2.0*sqrt(PHI + VP + 4.0*Vt));

            // Mobility equation
            double beta = KP* (W/L)* (1.0/(1.0 + THETA* VP));   

            // forward and reverse currents
            double xf = (VP-VS)/Vt;
            double iff = (log(1.0+exp( xf/2.0)))*(log(1.0+exp( xf/2.0)));
            double xr = (VP-VD)/Vt;
            double ir = (log(1.0+exp( xr/2.0)))*(log(1.0+exp( xr/2.0)));

            // Specific current
            double Ispec = 2 * n * beta * Vt * Vt;

            // Drain current
            double Id =0;

// VG>=VT0 进行Idot计算 
            Idot = set_final(VD,VG,Ipeak,Cd,Ctot,Beta);

        




            
    Id = Ispec * (iff - ir) + Idot;
    return Id;


}


double MNS(double a,double b, double c)
{
double d;
d = 0.5*(a+b-sqrt((a-b)*(a-b)+c));
return d;
}

double MXS(double a,double b, double c)
{
double d;
d = 0.5*(a+b+sqrt((a-b)*(a-b)+c));
return d;
}

};









#endif //__FUNC_H_
















