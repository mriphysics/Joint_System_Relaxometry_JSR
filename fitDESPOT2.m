 function [E_DESPOT2] = fitDESPOT2(noisyData,E_DESPOT1,rf_phase_incr,TRssfp,TE,alpha_ssfp,rf_trms,x0)
 noisyData  = double(noisyData);
 S          = @(E) (abs(simulateJSRvoxel(E(1),E_DESPOT1,E(2),E(3),...
     rf_phase_incr,0.001,TRssfp,TE,1,alpha_ssfp,rf_trms,true)));
 cf         = @(E) [S(E) - noisyData(:)];

 options                = optimset('lsqnonlin');
 options.FinDiffType    = 'forward';
 options.TolFun         = 1e-15;
 options.TolX           = 1e-15;
 options.MaxFunEvals    = 1e6;
 options.MaxIter        = 1e3;
 options.Display        = 'off';

 lb                     = [0 0 -inf];
 ub                     = [1e6 3000 inf];

 try
 [E_DESPOT2,~,~]        = lsqnonlin(cf,double(x0(1,:)),lb,ub,options);
 catch ME
     E_DESPOT2 = 0;
 end

 end