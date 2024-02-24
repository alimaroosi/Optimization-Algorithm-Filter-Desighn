%  Traning Feed-forward Neural Networks using Grey Wolf Optimizer   %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili,How effective is the Grey Wolf         %
%               optimizer in training multi-layer perceptrons       %
%              Applied Intelligece, in press, 2015,                 %
%               http://dx.doi.org/10.1007/s10489-014-0645-7         %
%                                                                   %

% Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
global Flag_Norm
global L_ave With_LM
% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
% % display('please comment line 34 to 36 of GWO just for test')
% % load ps3_PI1e5PF0.0157PT0.0117.mat
% % Positions(1,:)=hh;
Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter
% % % load Alpha_pos_PI
% % % Positions(1,:)= Alpha_pos;
% Main loop
while l<Max_iter
     LenDiv2= size(Positions,2)/2;
    for i=1:size(Positions,1)
        
        % Calculate objective function for each search agent
        [fitness x_Temp]=fobj(Positions(i,:));
        Positions(i,:)=x_Temp;
        if(Flag_Norm==1)
            Positions(i,:)=Positions(i,:)/sqrt(sum(Positions(i,:).^2));
        end
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness<Beta_score
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
        
        
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
            
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
            
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
            
            r1=rand();
            r2=rand();
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
                %%%%%%%%%%%Start movmean
        fitness=fitness1(Positions(i,:));
        x_Temp=Positions(i,:);
        X_movmean=[movmean(x_Temp(1:LenDiv2),L_ave) movmean(x_Temp(LenDiv2+1:end),L_ave)];
        [fitness_movmean X_movmean]=fitness1(X_movmean);
        if(fitness_movmean<fitness)
            fitness= fitness_movmean;
            Positions(i,:)=X_movmean;
        end
        %%%%%%%%%%%End movmean
        %%%%%%%%%LOcal
        T=l+1;
    end
    %%%%%%%%%LOcal 2
    T=l+1;
    if(With_LM && mod(T,floor(Max_iter/10))==0  && T>(3/4)*Max_iter)
        
        Tmax=Max_iter;
        options=optimoptions(@lsqnonlin);
        options.Algorithm='levenberg-marquardt';
        options.MaxFunctionEvaluations=2000;
        options.MaxIterations=1000;
        options.TolFun=1e-30;
        options.TolPCG=1e-30;
        options.TolX=1e-30;
        options.InitDamping=1e-20;
        %             if(T>2*Tmax/3)
        kk_Max=floor(30*(1-T/Tmax))+1;
        options.MaxFunctionEvaluations=T*floor(1000/kk_Max);
        options.MaxIterations=100000;
       
       
        
        x0=Beta_pos;
        [x_Temp] = lsqnonlin(fobj,x0,[],[],options);%,lb,ub);
        [fitness x_Temp]=fobj(x_Temp);
%         if(Flag_Norm==1)
%             x_Temp=x_Temp/sqrt(2*sum(x_Temp.^2)-x_Temp(1)^2);
%         end
        %%%%%%%%%%%Start movmean
        X_movmean=[movmean(x_Temp(1:LenDiv2),L_ave) movmean(x_Temp(LenDiv2+1:end),L_ave)];
        [fitness_movmean X_movmean]=fitness1(X_movmean);
        if(fitness_movmean<fitness)
            fitness= fitness_movmean;
            x_Temp=X_movmean;
        end
        %%%%%%%%%%%End movmean
        if( fitness<Beta_score)
            Beta_score=fitness;
            Beta_pos=x_Temp;
        end
        %%%
                x0=Alpha_pos;
        [x_Temp] = lsqnonlin(fobj,x0,[],[],options);%,lb,ub);
        [fitness x_Temp]=fobj(x_Temp);
%         if(Flag_Norm==1)
%             x_Temp=x_Temp/sqrt(2*sum(x_Temp.^2)-x_Temp(1)^2)
%         end

        %%%%%%%%%%%Start movmean
        X_movmean=[movmean(x_Temp(1:LenDiv2),L_ave) movmean(x_Temp(LenDiv2+1:end),L_ave)];
        [fitness_movmean X_movmean]=fitness1(X_movmean);
        if(fitness_movmean<fitness)
            fitness= fitness_movmean;
            x_Temp=X_movmean;
        end
        %%%%%%%%%%%End movmean

        if( fitness<Alpha_score)
            Alpha_score=fitness;
            Alpha_pos=x_Temp;
        end
        %%%
                x0=Delta_pos;
        [x_Temp] = lsqnonlin(fobj,x0,[],[],options);%,lb,ub);
        [fitness x_Temp]=fobj(x_Temp);
%         if(Flag_Norm==1)
%             x_Temp=x_Temp/sqrt(2*sum(x_Temp.^2)-x_Temp(1)^2)
%         end
        %%%%%%%%%%%Start movmean
        X_movmean=[movmean(x_Temp(1:LenDiv2),L_ave) movmean(x_Temp(LenDiv2+1:end),L_ave)];
        [fitness_movmean X_movmean]=fitness1(X_movmean);
        if(fitness_movmean<fitness)
            fitness= fitness_movmean;
            x_Temp=X_movmean;
        end
        %%%%%%%%%%%End movmean
        if( fitness<Delta_score)
            Delta_score=fitness;
            Delta_pos=x_Temp;
        end
        %%%
    end
    %%%%%%End Local2
    
    l=l+1;
    Convergence_curve(l)=Alpha_score;
    
    if(mod(l,50)==0)
        xout=Alpha_pos(1:length(Alpha_pos)/2);
        [PT,PF]=powerOOB(xout,64);
        [PI,~]=SMTINTR(xout,64);
        display (['Iteration = ' num2str(l) ' from ' num2str(Max_iter) ' PI= ' num2str(PI) 'PF= ' num2str(PF) ' PT= ' num2str(PT)]) 
    end
    
end



