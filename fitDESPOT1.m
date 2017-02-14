function [E_DESPOT1] = fitDESPOT1(noisyData,TRspgr,alpha_spgr)
%USAGE:
%    [E_DESPOT1] = fitDESPOT1(noisyData, TRspgr, alpha_spgr)
%
%    Estimate M_0 and T_1 via DESPOT1 method, data is linearized
%        by normalizing Ernst formula by sin(alpha_spgr)
%
%    noisyData  -> Measured SPGR data as is
%    TRspgr     -> SPGR repetition time
%    alpha_spgr -> excitation flip angle in radians
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk

Y = noisyData(:)./sin(alpha_spgr(:));
X = noisyData(:)./tan(alpha_spgr(:));
p = polyfit(X,Y,1);
m = p(1);
E_DESPOT1(2) = -TRspgr./log(m);
E_DESPOT1(1) = p(2)./(1-m);

end
