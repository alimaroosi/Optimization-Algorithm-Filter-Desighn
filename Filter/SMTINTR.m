function [PI,SIR]=SMTINTR(h,L)
%L=64;
%PI<10^(-4)
%Lh=257
h=[h(end:-1:2) h];
Npac=20; %number of packets data
h=h/sqrt(sum(abs(h.^2)));
Lh=length(h);
LE2=(Lh-1)/(L);
%Filter polyphase decomposition
for l=1:L
    hm=h(l:L:Lh);
    Lhm=length(hm);
    E(l,1:Lhm)=hm; %E_(l-1): (l-1)th polyphase decomposition filter of h
end
 


DATA1=zeros(L,Npac);
DATA1(L/2,Npac/2)=1;
%DATA1(2,1)=1;
DATA2=zeros(L,Npac);
% DATA2(6,3)=1;


l=1:L;
phi1=conj((1i).^(l-1))';
Phi1=repmat(phi1,1,Npac);
DATAp1=DATA1.*Phi1; %Product lth data to j^(l-1)

phi2=conj((1i).^(l))';
Phi2=repmat(phi2,1,Npac);
DATAp2=DATA2.*Phi2; %Product lth data to j^(l)

%Polyphase implementation of trasnmitter
x1=PSFB(DATAp1,L,E); %Pass data on n*T from polyphase synthesis filter bank 
x2=PSFB(DATAp2,L,E); %Pass data on n*T+T/2 from polyphase synthesis filter bank

x1h=[x1 zeros(1,L/2)];
x2h=[zeros(1,L/2) x2]; %L/2 time shift data on n*T+T/2 

x=x1h+x2h; %Transmitting data
Lx=length(x); %%Transmitting data length

y=x;

%Polyphase implementation of Receiver
YP1=PAFB(y,L,E); %Pass data from polyphase analysis filter bank to achieve data on n*T  
[LL,LYP1]=size(YP1);

y2=y(L/2+1:end); %-L/2 time shift data to achieve data on n*T+T/2
YP2=PAFB(y2,L,E); %Pass data from polyphase analysis filter bank to achieve data on n*T+T/2 
[LL,LYP2]=size(YP2);

l=1:L;
Nphi1=conj((-1i).^(l-1))';
NPhi1=repmat(Nphi1,1,LYP1);
YPP1=YP1.*NPhi1; %Product lth data to (-j)^(l-1)
YPP11=YPP1(:,LE2+1:LE2+1+Npac-1);

Nphi2=conj((-1i).^(l))';
NPhi2=repmat(Nphi2,1,LYP2);
YPP2=YP2.*NPhi2; %Product lth data to (-j)^(l)
YPP22=YPP2(:,LE2+1:LE2+1+Npac-1);


%Detection
YD1=(real(YPP11));
YD2=(real(YPP22));

PI=sum(sum(abs(YD1.^2)))+sum(sum(abs(YD2.^2)))-abs(YD1(L/2,Npac/2))^2;
SIR=abs(YD1(L/2,Npac/2))^2/PI;


