function [X_Best]=Run_PS(Alpha_pos, IterMax_PS)
warning('off','all')
addpath(genpath(cd))
% % % load Alpha_pos_pf016PI9e4PT010  %%% name of mat file resulted from GWO
global LFilter data
LFilter=length(Alpha_pos);
nvar=LFilter;

Abs_HB=Alpha_pos;
lb=Abs_HB-0.1*abs(Abs_HB);%movmin(Abs_HB,10)-0.01*abs(movmin(Abs_HB,10));%
ub= Abs_HB+0.1*abs(Abs_HB);%movmax(Abs_HB,10)+0.01*abs(movmax(Abs_HB,10));%
x0=Alpha_pos;
Tmax_PS=100;
data.lb=lb;
data.ub=ub;
X_Best=x0;

for i=1:IterMax_PS
    T=Tmax_PS/2;
    sol.x= X_Best;
    sol = fitness(sol,data);
    Fitness_Best=sol.fit;
    %%%%%
    [X_Best_Temp Fitness_Best_Temp]=PatternSearch( X_Best,Fitness_Best,nvar,data,  T, Tmax_PS);
    x=X_Best_Temp;
    
    
    
    %%%%%
    [fit_x x]=fitness1(x);
    [fit_X_Best X_Best]=fitness1(X_Best);
    
    
    if(fit_x< fit_X_Best)
        X_Best=x;
    end
    [PT,PF]=powerOOB(X_Best(1:LFilter/2),64);
    [PI,SIR]=SMTINTR(X_Best(1:LFilter/2),64);
   display (['external Iteration = ' num2str(i)  ' PI= ' num2str(PI) 'PF= ' num2str(PF) ' PT= ' num2str(PT) 'SIR= ' num2str(SIR)]) 
    save X_outBest X_Best;
end
