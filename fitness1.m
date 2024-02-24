function [O xout]=fitness1(x)
global data
global LFilter
global With_Freq 
global With_Time 

% % global shift_x
% % O=sum((x-shift_x).^2);
sol1.x=x(1:LFilter/2);
data_temp.lb=data.lb(1:LFilter/2);
data_temp.ub=data.ub(1:LFilter/2);
[sol1,O]=fitness(sol1,data_temp);
Hx=x(LFilter/2+1:end);
Hx=[Hx Hx(end:-1:2)];
hx=ifft(Hx);
sol2.x=hx(1:LFilter/2);
[sol2,O]=fitness(sol2,data_temp);
if(With_Freq ==0)
sol2.fit=1e5; %%% remove frequency part if enable 
end
if(With_Time==0)
   sol1.fit =1e5;
end
if(sol1.fit<=sol2.fit)
O=sol1.fit;
  x_Temp=sol1.x;
  x_Temp=x_Temp/sqrt(2*sum(x_Temp.^2)-x_Temp(1)^2);
    H=ifft([x_Temp x_Temp(end:-1:2)]);
else
  O=sol2.fit;
  x_Temp=sol2.x;
  x_Temp=x_Temp/sqrt(2*sum(x_Temp.^2)-x_Temp(1)^2);
  H=ifft([x_Temp x_Temp(end:-1:2)]);
end
xout=[x_Temp H(1:LFilter/2)];