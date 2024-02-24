function [PT,PF, H]=powerOOB(h,L)
%L=64
%PF minimize
%PF approximately is 0.01
global PF_Weight

h=[h(end:-1:2) h];
h=h/sqrt(sum(abs(h.^2)));
Lh=length(h);
LE=((Lh-1))/L;

T=sqrt(2); %Space between two symboles in time domain (FMT)
Ts=T/L;
Qh=1; %Oversampling rate;
TsQ=Ts/Qh; %Oversampling time duration
F=2/T; %Space between carriers in FMT

N0h=fix(T/(2*TsQ));
PDFT=2*LE*N0h;
Fs=1/(PDFT+1);
K0h=fix(F*TsQ/(2*Fs));


Lh2=(Lh-1)/2;
nn=-Lh2:Lh2;
W=exp(-j*2*pi*nn'*nn/Lh);
H=real(W*conj(h'));
H=H/sqrt(sum(abs(H.^2)));

if (PF_Weight==1) 
 PF=15*sum(abs(H(1:Lh2+1-3*K0h).^2))+4*sum(abs(H(Lh2+1-3*K0h+1:Lh2+1-2*K0h).^2))+sum(abs(H(Lh2+1-2*K0h+1:Lh2+1-K0h).^2))+sum(abs(H(Lh2+1+K0h:Lh2+1+2*K0h-1).^2))+4*sum(abs(H(Lh2+1+2*K0h:Lh2+1+3*K0h-1).^2))+15*sum(abs(H(Lh2+1+3*K0h:end).^2));
% % PF=1e10*sum(abs(H(1:Lh2+1-3*K0h).^2))+4*sum(abs(H(Lh2+1-3*K0h+1:Lh2+1-2*K0h).^2))+sum(abs(H(Lh2+1-2*K0h+1:Lh2+1-K0h).^2))+sum(abs(H(Lh2+1+K0h:Lh2+1+2*K0h-1).^2))+4*sum(abs(H(Lh2+1+2*K0h:Lh2+1+3*K0h-1).^2))+1e10*sum(abs(H(Lh2+1+3*K0h:end).^2));
% % tt=0;
elseif (PF_Weight==2)
	PF=4*sum(abs(H(1:Lh2+1-2*K0h).^2))+sum(abs(H(Lh2+1-2*K0h+1:Lh2+1-K0h).^2))+sum(abs(H(Lh2+1+K0h:Lh2+1+2*K0h-1).^2))+4*sum(abs(H(Lh2+1+2*K0h:end).^2));
else 
PF=sum(abs(H(1:Lh2+1-K0h).^2))+sum(abs(H(Lh2+1+K0h:end).^2));
end

PT=sum(abs(h(1:Lh2+1-N0h).^2))+sum(abs(h(Lh2+1+N0h:end).^2));
tt=0;