clc
clear

Vgs=0.44:0.0002:0.55;

Vds=0.006;
% Vgs= 0.4932;

temp=4.2;
qc = -1.60217733e-19;
kB=1.380658e-23;
m0=9.10956e-31;
weff=55e-9;
h=6.62607015e-34;

meff=0.19*m0;


Vpeak=[0.4704 0.4932 0.5136];
Ipeak=[4.2076E-09 3.42005e-08 4.53746e-08];
beta=[77.14772727];

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
s2=(0.4866-0.4924)/(0.0117-0);
s1=(0.4866-0.4704)/(0.0117-0);


%% caculate parameter such as Cg
Cg=abs(qc/((Vpeak(2)-Vpeak(1))-(s2-s1)/(s1*s2)*(delta_E(1)/qc)));
Cd=abs(Cg/s2);
Ctot=Cg*abs((s2-s1)/(s1*s2));
%% caculate energy of tunneling


It=0;

for i=1:3
    if i==1
        Qi=0;
    else
        Qi=0;
    end
    
    delta_FS_i_sub = qc/Ctot*(Cd*Vds + Cg*(Vgs-Vpeak(i))+Qi) + abs(delta_E(i));
    delta_FS_isub1_add = - delta_FS_i_sub;
    delta_FD_i_sub = qc/Ctot*(-(Ctot-Cd)*Vds+Cg*(Vgs-Vpeak(i))+Qi)+abs(delta_E(i));
    delta_FD_isub1_add = - delta_FD_i_sub;

    I_vgs = Ipeak(i)*exp(beta*(Vgs-Vpeak(i)));



    It = It +  I_vgs.*(1-exp(qc*Vds/(kB*temp)))./(1 +  exp(delta_FS_i_sub/(kB*temp))...
                                           +exp(-delta_FD_i_sub/(kB*temp))...
                                           +exp(qc*Vds/(kB*temp)) );
end
                                       
plot(Vgs,(It))



