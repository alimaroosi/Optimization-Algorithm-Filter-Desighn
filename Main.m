clc
clear
close all
warning('off','all')
format shortG
global FuncType_Global
global shift_x
global Flag_Norm
FuncType_Global=1; %1 run usual 2
global LFilter
global LB_Glob UB_Glob  data
global Thr_PT Thr_PF Thr_PI
global L_ave With_Freq With_Time With_LM With_PS
global PF_Weight

addpath(genpath(cd))

%% Parameters Definiterion
PF_Weight=1; %% 0 =no weight for lobe   1 =type1     2=Type2 
Flag_Norm=1; %%% normalize energy to 1
With_Smooth=1; %%% if 1  smoothing filte else no
With_Freq=1; %%% if 1 frequency domain also considered if 0 not
With_Time=1; %%% if 1 time domain also considered if 0 not
With_LM=1; %%% if 1 LevenbergMarquate applied if 0 not
With_PS=1; %%% if 1  pattern search applied else no
HalfFilterPlusone=257; %% 129 or 257, Half of filter length +1 should be 129 for 257 and 257 for 513
Thr_PT=0.0582; %% e.g. 0.058 prefered PT   For iota 0.0117
Thr_PF=1e-6;%0.005; %% e.g. 1e-7 Prefered PF     For iota  0.022
Thr_PI=3e-8; %% e.g. 1e-7 Prefered PI      For iota  1e-4
IterMax_PS=10;  %%% atleast 20 Maximum Number of iterations for PS default should be 20
maxiter_GWO=1000; % at least 1000 "Maximum Number of iterations for GWO" 
SearchAgents_no=30; %%%default should be 30
if(With_Smooth)
  L_ave=3; %%% should be 3 to 5 (Length of moving average)
else
  L_ave=1; %%% 1 means no moving average
end

if(~With_Freq & ~With_Time)
    display('error:atleast one of With_Freq or With_Time should be 1'); 
return
end
LFilter=HalfFilterPlusone*2; %% consider both time nad frequency domain
nvar=LFilter; %%% half of a solution dimension for time samples others for frequencies
% load hSM257.mat
% load 128_Pf012PI5e4PT010.mat
% if(h(1)>0.5)
%     h=[h(end:-1:2) h];
% end
% Abs_HB=(h(LFilter:end));
% lb=movmin(Abs_HB,10)-0.2*abs(movmin(Abs_HB,10));%Abs_HB-0.4*abs(Abs_HB);
% ub= movmax(Abs_HB,10)+0.2*abs(movmax(Abs_HB,10));%Abs_HB+0.4*abs(Abs_HB);%

lb=[-0.2*ones(1, nvar/2) -sqrt(LFilter)*ones(1, nvar/2)];
ub=[1*ones(1, nvar/2) sqrt(LFilter)*ones(1, nvar/2)];
data.lb=lb;
data.ub=ub;
LB_Glob=lb;
UB_Glob=ub;
Max_iter=maxiter_GWO;
dim=nvar;
fobj=@fitness1;
[Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj);
% load Alpha_pos_hsm_pf6e8PI1e7PT058.mat
h=Alpha_pos(1:LFilter/2);
h_After_GWO=h;
save h_After_GWO_Data h_After_GWO
[PT,PF]=powerOOB(h,64);
[PI,SIR]=SMTINTR(h,64);
display (['final after GWO h= ' num2str(h)]);
display (['Final Result after GWO '  ' PI= ' num2str(PI) ' PF= ' num2str(PF) ' PT= ' num2str(PT) 'SIR= ' num2str(SIR)])

if(With_PS)
[X_Best]=Run_PS(Alpha_pos, IterMax_PS);
[PT,PF]=powerOOB(X_Best(1:LFilter/2),64);
[PI,SIR]=SMTINTR(X_Best(1:LFilter/2),64);

h=X_Best(1:LFilter/2);
h_After_PS=h;
save h_After_PS_Data h

display (['final after PS h= ' num2str(h)]);
display (['Final Result '  ' PI= ' num2str(PI) ' PF= ' num2str(PF) ' PT= ' num2str(PT) ' SIR= ' num2str(SIR)])

plot(h)
else
    display('Note: PS not applied')
end
cc=0;