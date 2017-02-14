function [E_DESPOT2] = fitDESPOT2(noisyData,E_T1,rf_phase_incr,TRssfp,TE,alpha_ssfp,rf_trms,x0)
% USAGE:
%    [E_DESPOT2] = fitDESPOT2(noisyData,E_T1,rf_phase_incr,TRssfp,TE,alpha_ssfp,rf_trms,x0)
%
%    Estimate M_0, T_2 and B_0 via DESPOT2-FM method, note that here we apply levenberg-marquardt
%        algorithm instead of stocastic region contraction that was first proposed.
%    Assumes T_1 was previously estimated in a seperate measurement.
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk

noisyData            = double(noisyData);

options              = optimoptions('lsqnonlin');
options.Algorithm    = 'levenberg-marquardt';
options.TolFun       = 1e-15;
options.TolX         = 1e-15;
options.MaxFunEvals  = 1e6;
options.MaxIter      = 1e3;
options.Display      = 'off';
options.Jacobian     = 'on';
lb                   = [];
ub                   = [];
%lb = [0 0 -inf];
%ub = [1e6 3000 inf];

try
    [E_DESPOT2,~,~]  = lsqnonlin(@(E) despot2fm_cost_function(noisyData(:),E,),double(x0(1,:)),lb,ub,options);
catch ME
    E_DESPOT2        = 0;
end
end
