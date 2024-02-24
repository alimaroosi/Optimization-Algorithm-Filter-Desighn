function x=PSFB(DATAp,L,E)

%Polyphase Synthesis Filter Bank
%DATAp: Input data
%L: Number of subcarriers
%E: Transmitter polyphase decomposition filters 

%Step 1: IDFT of input Data
data=L*ifft(DATAp,L,1);

%Step 2: pass data from synthesis polyphase filters
for l=1:L
    xP(l,:)=conv(data(l,:),E(l,:));
end

%Step 3: Parallel to serial data
[LL,Ld]=size(xP);
x=[];
for l=1:Ld
    x=[x conj(xP(:,l)')];
end