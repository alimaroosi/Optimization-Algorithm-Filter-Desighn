function [sol,O]=fitness(sol,data)
global Thr_PT Thr_PF Thr_PI
x=sol.x;

lb=data.lb;
ub=data.ub;
x=CB(x,lb,ub);
sol.x=x;
L= 64;
h=x;
[PT,PF]=powerOOB(h,L);  %%PF minimize 0.01
[PI,~]=SMTINTR(h,L);  %%%PI<1e-4
Level=1e-4;

%%%%%For hsm 257
w_PF=Thr_PI/(Thr_PI+Thr_PF);
w_PI=1-w_PF;
Beta=10;
if(PT>Thr_PT)
% %     fit=5+w_PF*PF+w_PI*PI+Beta*PT;
    fit=5+1*PF+1*PI+Beta*PT;
elseif(PT<(Thr_PT-0.1*Thr_PT))
% %     fit=3+w_PF*PF+w_PI*PI+ abs(PT-Thr_PT);
    fit=3+1*PF+1*PI+ abs(PT-Thr_PT);

elseif(PF < Thr_PF & PI <Thr_PI)
    fit=w_PF*PF+w_PI*PI;
else
    if(PF > 5*max(Thr_PF,Thr_PI) || PI>5*max(Thr_PF,Thr_PI)) %%%prevent from bad pattern of shape
      fit=2.5+1*PF+1*PI;
    elseif(PF < 0.9*Thr_PF)
     fit=1.5+w_PF*PF+w_PI*PI+10*abs(PF-Thr_PF);
    elseif(PI < 0.9*Thr_PI)
     fit=1.5+w_PF*PF+w_PI*PI+10*abs(PI-Thr_PI);
     else
      fit=w_PF*PF+w_PI*PI;  
    end
end
%%%%%End hsm 257
O=fit;
sol.fit=fit;

end





function x=CB(x,lb,ub)

x=max(x,lb);
x=min(x,ub);

end
