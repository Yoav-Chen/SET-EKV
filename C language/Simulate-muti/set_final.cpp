#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>

#include "set_final.h"

#define eps 1e-40

double set_final(double Vds,double Vgs,double *Ipeak,double *Cd,double *Ctot,double *beta){


    double temp = 4.2;
    double qc = -1.60217733e-19;
    double kB=1.380658e-23;

    double m0=9.10956e-31;
    double weff = 55e-9;
    double h= 6.62607015e-34;
    double meff=0.19*m0;

    double delta_E[6] = {0e-21, 0.3145e-21, 0e-21, 0.5241e-21, 0e-21, 0.7337e-21};
    double Cg[3] = {0,0,0};
/*     double Cd[3] = {0,0,0};
    double Ctot[3] = {0,0,0};

    double Ipeak[3] = {0,0,0};
    double beta[3] = {0,0,0}; */

    double delta_FS_i_sub = 0;
    double delta_FS_isub1_add = 0;
    double delta_FD_i_sub = 0;
    double delta_FD_isub1_add = 0;
        
    double I_vgs = 0;
    double Idot = 0;

    double Vpeak[3] = {0.4704, 0.4924, 0.5136};


    int i;

    for (i=2;i>=0;i--){
        if (i!=0)
        {
            Cg[i]=abs(qc/(Vpeak[i]-Vpeak[i-1]));
        }    
        else 
        {
            Cg[0]=Cg[1];
        }
    }



    for (i=0;i<3;i++){
        
        delta_FS_i_sub = qc/Ctot[i]*(Cd[i]*Vds + Cg[i]*(Vgs-Vpeak[i])) + abs(delta_E[i]);
        delta_FS_isub1_add = - delta_FS_i_sub;
        delta_FD_i_sub = qc/Ctot[i]*(-(Ctot[i]-Cd[i])*Vds+Cg[i]*(Vgs-Vpeak[i]))+abs(delta_E[i]);
        delta_FD_isub1_add = - delta_FD_i_sub;

        I_vgs = Ipeak[i]*exp(beta[i]*(Vgs-Vpeak[i]));





        Idot = Idot + I_vgs * (1-exp(qc*Vds/(kB*temp)) )/ (1 + exp(delta_FS_i_sub/(kB*temp)) \
                                                             + exp(-delta_FD_i_sub/(kB*temp)) \
                                                             + exp(qc*Vds/(kB*temp))    );


    }
    

    return Idot;
}

























