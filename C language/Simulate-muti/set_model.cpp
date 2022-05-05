#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>

#include "set_model.h"

#define eps 1e-40

double set_model(double Vd,double Vs,double Vg,double Vb,double Vpeak,double temp,int number,double Cd,double Cs,double Cg,double Cb,double Rd,double Rs,double C0,double Q0,double k){
    double eV = 1.60217733e-19;
    double kB=1.380658e-23;
    double Idot=0;
    double ref_vg=0;
    double distance=0;
    double Csum = Cd+Cs+Cg+Cb;
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

    int i;

    int i_max=0;

    i_max = round( -(Q0 + Cd*Vd + Cs*Vs + Cg*(Vg-Vpeak) + Cb*Vb)/eV + Csum/eV * (Rd*Vs + Rs*Vd)/(Rd + Rs) + 0.5 );
    


    P[(number-1)/2]=1;

    
    for (i=0;i<number;i++){
        i_state[i] = i -  ((number-1)/2 - i_max); 
        Vi[i] = (i_state[i]*eV + Q0 + Cd*Vd + Cs*Vs + Cg*(Vg-Vpeak) + Cb*Vb) / Csum;
        d_Fdi[i] =  -eV*(Vd-Vi[i])+ pow(eV,2) / (2*Csum);
        d_Fid[i] =  -eV*(Vi[i]-Vd)+ pow(eV,2) / (2*Csum);
        d_Fsi[i] =  -eV*(Vs-Vi[i])+ pow(eV,2) / (2*Csum);
        d_Fis[i] =  -eV*(Vi[i]-Vs)+ pow(eV,2) / (2*Csum);
        Gamma_DR[i] =  -d_Fdi[i] / (pow(eV,2) * Rd * (1-exp(d_Fdi[i]/(kB*temp))));
        Gamma_DL[i] =  -d_Fid[i] / (pow(eV,2) * Rd * (1-exp(d_Fid[i]/(kB*temp))));
        Gamma_SL[i] =  -d_Fsi[i] / (pow(eV,2) * Rs * (1-exp(d_Fsi[i]/(kB*temp))));
        Gamma_SR[i] =  -d_Fis[i] / (pow(eV,2) * Rs * (1-exp(d_Fis[i]/(kB*temp))));
    }

    //  for (i=0;i<number;i++){
    //      printf("%0.2e %0.2e %0.2e %0.2e \n", Gamma_DR[i],Gamma_DL[i], Gamma_SL[i],Gamma_SR[i] );
    //  }

    for (i=(number-1)/2+1;i<number;i++){
        P[i]=P[i-1]*( Gamma_DR[i-1] + Gamma_SL[i-1] ) / (Gamma_DL[i] + Gamma_SR[i]);
    }
    
    for (i=(number-1)/2-1;i>=0;i--){
        P[i]=P[i+1]*( Gamma_DL[i+1] + Gamma_SR[i+1] ) / (Gamma_DR[i] + Gamma_SL[i]);
    }

    
   

    for (i=0;i<number;i++){
        P_sum = P_sum + P[i];
    }

    for (i=0;i<number;i++){
        P[i] = P[i]/P_sum;
    }    



    for (i=0;i<number;i++){
        Idot = Idot+eV*P[i]*(Gamma_DR[i]-Gamma_DL[i]);
    }

/*     for (i=0;i<number;i++){
        Idot = Idot+eV*P[i]*(Gamma_DR[i]-Gamma_DL[i])*amp[i_max];
    } */
/*     if(-Cg/Cd*(Vg-(Vpeak+0.5*eV/Cg))-Vd > eps && Vd>0)
    {
        ref_vg = Vd*(-Cd/Cg)+(Vpeak+1/2*eV/Cg);
        distance = ref_vg-Vg;
        Idot = (1-tanh(k*distance))*set_model(Vd,Vs,ref_vg,Vb,Vpeak,temp,number,Cd,Cs,Cg,Cb,Rd,Rs,C0,Q0,k);
    }
        
    else if (Cg/(Cg+Cs)*(Vg-(Vpeak+0.5*eV/Cg))-Vd < -eps && Vd<0)
    {
        ref_vg = Vd*((Cg+Cs)/Cg) + (Vpeak+1/2*eV/Cg);
        distance = ref_vg-Vg;
        Idot = (1-tanh(k*distance))*set_model(Vd,Vs,ref_vg,Vb,Vpeak,temp,number,Cd,Cs,Cg,Cb,Rd,Rs,C0,Q0,k);
    }
     */

    return Idot;

}