clc
clear

% Vgs=0.44:0.0002:0.525;

Vds=0.002;
Vgs= 0.4704;

temp=4.2;
qc = -1.60217733e-19;
kB=1.380658e-23;
m0=9.10956e-31;
weff=55e-9;


h=6.62607015e-34;

meff=0.19*m0;


Vpeak=[0.4704 0.4924 0.5136];
% Ipeak=[4.2076E-09 3.42005e-08 4.53746e-08];  %original
Ipeak=[4.30021e-09 1.56846e-08 2.24052e-08];
beta=[75.6174  60.985 -50.3666];

delta_E=zeros(1,6);
for i=1:6
    if mod(i,2)==1  %奇数
        delta_E(i)=0;
    else            %偶数
        delta_E(i)=(i+1)*h^2/(8*meff*weff^2);
    end
end

i=1;

%% slope of blockade
s2=(0.487-0.4924)/(0.0112-0);
s2=2.5*s2
s1=(0.487-0.4704)/(0.0112-0);
s1=0.45*s1

%% caculate parameter such as Cg
% Cg=abs(qc/((Vpeak(2)-Vpeak(1))-(s2-s1)/(s1*s2)*(delta_E(1)/qc)));
% Cd=abs(Cg/s2);
% % Cd=10e-18;
% Ctot=Cg*abs((s2-s1)/(s1*s2));


for i=3:-1:1
    if( i~=1 )
        Cg(i)=abs(qc/(Vpeak(i)-Vpeak(i-1)));
    elseif(i==1)
        Cg(1)=Cg(2);
    end
    
end

Cd = [6.85626e-18 7.0545e-18 7.54572e-18];

Ctot=[1.74584e-17 1.45197e-17 1.89818e-17];


%% caculate energy of tunneling


It=0;

for i=1:1
    if i==1
        Qi=0;
    else
        Qi=0;
    end
    
    delta_FS_i_sub = qc./Ctot(i).*(Cd(i)*Vds + Cg(i).*(Vgs-Vpeak(i))+Qi) + abs(delta_E(i));
    delta_FS_isub1_add = - delta_FS_i_sub;
    delta_FD_i_sub = qc/Ctot(i).*(-(Ctot(i)-Cd(i))*Vds+ Cg(i).*(Vgs-Vpeak(i))+Qi)+abs(delta_E(i));
    delta_FD_isub1_add = - delta_FD_i_sub;

    I_vgs = Ipeak(i).*exp(beta(i).*(Vgs-Vpeak(i)));



    It = It +  I_vgs.*(1-exp(qc*Vds/(kB*temp)))./(1 +  exp(delta_FS_i_sub/(kB*temp))...
                                           +exp(-delta_FD_i_sub/(kB*temp))...
                                           +exp(qc*Vds/(kB*temp)) );
                                       
    
end

mode =1;   % @1 original data @ 2 sub EKV
if mode == 1
    It = It + EKV_model(Vds,Vgs);
elseif mode == 2
    ;
end


def = 2;   % @1 log     @2 linear
if def ==1
   plot(Vgs,log(It+1E-10)/log(10))
elseif def == 2
    plot(Vgs,abs(It))
end





