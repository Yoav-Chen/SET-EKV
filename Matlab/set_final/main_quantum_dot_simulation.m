clc
clear

Vs=0;
Vb=0;

Vd=-0.006:0.0001:0.006;
Vg=0.44:0.0002:0.525;

Ids=[];

for vd = Vd
    Ids_vg=[];
    imax_vg=[];
    for vg = Vg
       [current] = function_quantum_dot_simulation(vd,vg)+EKV_model(vd,vg);
% [current] = EKV_model(vd,vg);
       Ids_vg=[Ids_vg,current];
    end
    Ids=[Ids;Ids_vg];

end

[Vgs,Vds]=meshgrid(Vg,Vd);

mesh(Vgs,Vds,(abs(Ids)+1E-13))

% mesh(Vgs,Vds,log(abs(Ids)+1E-10)/log(10))
% set(gcf, 'Renderer', 'opengl')
% set(gcf,'paperpositionmode','auto')
% print('-depsc','a.eps')
