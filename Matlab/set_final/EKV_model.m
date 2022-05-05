function [Ids]=EKV_model(VD,VG)
    
    W=0.3;
    L=0.055;
    
    Vt = 0.00010341;
    VTO = 0.483933;
    GAMMA = 0.00674447;
    PHI = 1.77018e+12;
    
    
%     KP = 5.36202e-05;
    KP = 4.5e-05;
    THETA = 6.44918e-08;
    
    VB=0;
    VS=0;
    VGprime = VG - VTO + PHI + GAMMA * sqrt(PHI);
    VP = VGprime - PHI - GAMMA * (sqrt(VGprime+(GAMMA/2.0)*(GAMMA/2.0))-(GAMMA/2.0));
    
    n = 1.0 + GAMMA./(2.0*sqrt(PHI + VP + 4.0*Vt));
    beta = KP* (W/L)* (1.0./(1.0 + THETA.* VP)); 
    xf = (VP-VS)/Vt;
    iff = (log(1.0+exp( xf/2.0))).*(log(1.0+exp( xf/2.0)));
    xr = (VP-VD)./Vt;
    ir = (log(1.0+exp( xr/2.0))).*(log(1.0+exp( xr/2.0)));
    Ispec = 2 .* n .* beta * Vt * Vt;
    Ids = Ispec .* (iff - ir);
    
end
