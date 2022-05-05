clc
clear

Vd=[0 0.001 0.002 0.003 0.004 0.005 0.006];
Vg=0.35:0.0002:0.52;



for j=1:length(Vd)
    Ids=[];
    for i= 1:length(Vg)
       Ids=[Ids,EKV_model(Vd(j),Vg(i))];
    end

plot(Vg,(Ids))

hold on
end

