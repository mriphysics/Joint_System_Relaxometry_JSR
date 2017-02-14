%% Function that estimates M0,T1,T2,B1 and B0 based on the JSR approach

%
% Rui Pedro A. G. Teixeira @ 31/05/2016
% Rui Pedro A. G. Teixeira @ 30/06/2016 - allow bounds to be set externaly

function [E_parameters, resNorm, res, exitFlag] = fit_jsr_model(noisySignal,...
        RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr,...
        x0,options,lb,ub)
    
    noisySignal     = double(noisySignal);
    
    %if user doesn't give a starting point we do an 'educated guess'
    if isempty(x0)
        x0 = [10*noisySignal(1) 0 900 50 0];
        
    end
       
    
    % if options are not given as input create option structure
    if isempty(options)
        
        options = optimset('lsqnonlin');
        options.FinDiffType = 'central';
        options.TolFun = 1e-4;
        options.TolX = 1e-4;
        options.MaxFunEvals     = 1e12;
        options.MaxIter = 1e12;
        options.Display = 'off';
        
    end
    %

    [E_parameters,resNorm,res,exitFlag] = lsqnonlin(@(E) jsr_cost_function(noisySignal(:),E,RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr),x0,lb,ub,options);

end