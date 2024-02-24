function YP=PAFB(y,L,E)

%Polyphase Analysis Filter Bank
%YP: Receiving data
%L: Number of subcarriers
%E: Transmitter polyphase decomposition filters 

Ly=length(y); %Receiving data length

%Step 1: Serial to parallel
ym=y(1:L:Ly);
Lym=length(ym);
rP(1,1:Lym)=ym;
for l=2:L
    ym=y(L-(l-2):L:Ly);
    Lym=length(ym);
    rP(l,2:Lym+1)=ym;
end
% Ld=floor(Ly/L);
% for l=1:Ld
%     ym=conj(y((l-1)*L+1:l*L)');
%     Lym=length(ym);
%     rP(:,l)=ym;
% end
%Step 2: Pass data from analysis polyphase filters
for l=1:L
    yP(l,:)=conv(rP(l,:),E(l,:));
end

%Step 3: IDFT
YP=L*ifft(yP,L,1);