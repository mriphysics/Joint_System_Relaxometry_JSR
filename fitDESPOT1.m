 function [E_DESPOT1] = fitDESPOT1(noisyData,TRspgr,alpha_spgr)

 Y = noisyData(:)./sind(alpha_spgr(:));
 X = noisyData(:)./tand(alpha_spgr(:));
 p = polyfit(X,Y,1);
 m = p(1);
 E_DESPOT1(2) = -TRspgr./log(m);
 E_DESPOT1(1) = p(2)./(1-m);
 end
