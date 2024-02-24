function  x_Temp=Filterd(x_Temp)
global h_rcs lpFilt
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.02, ...
%          'StopbandFrequency',0.05,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% fvtool(lpFilt)
% x_Temp=h_rcs(257:end);

% [b,a] = ellip(6,1,70,1/513);
%  dataOut = filter(b,a,[h_rcs h_rcs h_rcs]);
%     y=filter(lpFilt,[ x_Temp(end:-1:2) x_Temp x_Temp(end:-1:2) x_Temp]);
%     start_sel=length(x_Temp)+lpFilt.filtord/2;
%     x_Temp=y(start_sel:start_sel+256);
x_Temp=x_Temp;
end