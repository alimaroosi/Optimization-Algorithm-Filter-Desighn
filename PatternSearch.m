function [X_Best Fitness_Best]=PatternSearch( X_Best,Fitness_Best,Dimension,data,  T, Tmax)
global L_ave
Max_minus_Min=data.ub-data.lb;


%%%%%%%% New local search only best
LenDiv2=length(X_Best)/2;
q=0;
for Try=1:10
    
    ChosenDim= randperm(Dimension); %%% dimension chosen randomly
    for i =1:Dimension
        Patern_Dim= ChosenDim(i);
        PaternDel= Max_minus_Min(Patern_Dim)/4;
        
        
        while PaternDel> (1-T/Tmax)^2*Max_minus_Min(Patern_Dim)* .01  %0.0001
            X_Patern=X_Best;
            X_Patern(Patern_Dim)= X_Best(Patern_Dim)+ PaternDel;
            sol.x=X_Patern;
            [sol.fit sol.x]=fitness1(X_Patern);
            Fitness_Patern= sol.fit;
            %%%%%%%%%%%Start movmean
            X_Patern2=[movmean(X_Patern(1:LenDiv2),L_ave) movmean(X_Patern(LenDiv2+1:end),L_ave)];
            [sol2.fit sol2.x]=fitness1(X_Patern2);
            if(sol2.fit<sol.fit)
                Fitness_Patern= sol2.fit;
                X_Patern=sol2.x;
            end
            %%%%%%%%%%%End movmean
            
            if(Fitness_Patern<Fitness_Best)
                X_Best = X_Patern;
                Fitness_Best = Fitness_Patern;
            else
                X_Patern(Patern_Dim)= X_Best(Patern_Dim)- PaternDel;
                
                [sol.fit sol.x]=fitness1(X_Patern);
                Fitness_Patern= sol.fit;
                %%%%%%%%%%%Start movmean
                X_Patern2=[movmean(X_Patern(1:LenDiv2),L_ave) movmean(X_Patern(LenDiv2+1:end),L_ave)];
                [sol2.fit sol2.x]=fitness1(X_Patern2);
                if(sol2.fit<sol.fit)
                    Fitness_Patern= sol2.fit;
                    X_Patern=sol2.x;
                end
                %%%%%%%%%%%End movmean
                if (Fitness_Patern<Fitness_Best)
                    X_Best = X_Patern;
                    Fitness_Best = Fitness_Patern;
                else
                    PaternDel = PaternDel/2;
                end
                
                
            end
        end
            if(mod(q,50)==0)
                [PT,PF]=powerOOB(X_Best(1:LenDiv2),64);
               [PI,SIR]=SMTINTR(X_Best(1:LenDiv2),64);
               display (['Inter Iter PS= ' num2str(q) ' of ' num2str(Dimension*10)  ' PI= ' num2str(PI) 'PF= ' num2str(PF) ' PT= ' num2str(PT)]) 
                h=X_Best(1:LenDiv2);
                save h_After_PS_Data h
            end
            q=q+1;
    end
    
end

