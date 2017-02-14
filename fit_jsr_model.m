function [E_parameters, resNorm, res, exitFlag] = fit_jsr_model(noisySignal,...
        RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr,...
        x0,options,lb,ub)
% Function that estimates rM0, iM0, T1,T2, and B0 based on the JSR approach
%
%USAGE:
%     [E_parameters, resNorm, res, exitFlag] = fit_jsr_model(noisySignal,...
%     RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr,...
%     x0,options,lb,ub)
%
%    noisySignal -> expected as the concatenation of SPGR and rSSFP and iSSFP signals -  [SPGR1;...;SPGRn;rSSFP1;...;rSSFPn;iSSFP1;...;iSSFPn]
%    RF_dur       -> Duration of the rf pulse used to produce the excitation (ms)
%    T_E          -> Echo time applied, typically 0.5 T_R_SSFP (ms)
%    T_R_spgr     -> Repetition time of SPGR signal (ms)
%    T_R_ssfp     -> Repetition time of ssfp signal (ms)
%    a_spgr       -> Flip angles used in SPGR acquisition (rad)
%    a_ssfp       -> Flip angles used in SSFP acquisition (rad)
%    phi_incr     -> SSFP rf phase increment that allows to "move" characteristic bands
%    x0           -> Inital starting point for the search algorithm, if empty the default value of [10*noisySignal(1) 0 900 50 0] is used
%    options      -> matlab optimset structure:
%                                               options             = optimset('lsqnonlin');
%                                               options.TolFun      = 1e-15;
%                                               options.TolX        = 1e-15;
%                                               options.MaxFunEvals = 1e4;
%                                               options.MaxIter     = 1e4;
%                                               options.Display     = 'off';
%                                               options.Jacobian    = 'on';
%    lb           -> optional lower limit of the search space (necessary to be [] of LM algorithm is to be used)
%    ub           -> optional higher limit of the search space (necessary to be [] of LM algorithm is to be used)
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk
    
    noisySignal     = double(noisySignal); % lsqnonlin requires double precision data
    
    %if user doesn't give a starting point we do an 'educated guess'
    if isempty(x0)
        x0 = [10*noisySignal(1) 0 900 50 0];
        
    end
       
    
    % if options are not given as input create option structure
    if isempty(options)
        
        options             = optimset('lsqnonlin');
        options.TolFun      = 1e-15;
        options.TolX        = 1e-15;
        options.MaxFunEvals = 1e4;
        options.MaxIter     = 1e4;
        options.Display     = 'off';
        options.Jacobian    = 'on';
        
    end
    
    % perform the actual fit
    [E_parameters,resNorm,res,exitFlag] = lsqnonlin(@(E) jsr_cost_function(noisySignal(:),E,RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr),x0,lb,ub,options);

end
