function [model] = simulateJSRvoxel_complexCF(M0,T1,T2,b0,rf_phase_incr,TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms)
% Funtion that allows alocates the SPGR and bSSFP signals in a single array such that a joint estimation can be performed
% always outputs in an order such that:
%   model = [SPGR1; SPGR2; ..., SPGRNspgr;real(bSSFP1);real(bSSFP2);...;real(bSSFPNssfp);imag(bSSFP1);imag(bSSFP2);...;imag(bSSFPNssfp)];
%
% Usage:
%   [model] = simulateJSRvoxel_complexCF(M0,T1,T2,b0,rf_phase_incr,TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms)
%
Nspgr = length(alpha_spgr);

S = simulateJSRvoxel(M0,T1,T2,b0,rf_phase_incr,...
    TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms,false);

model = [abs(S(1:Nspgr));real(S(Nspgr+1:end));imag(S(Nspgr+1:end))];

end
