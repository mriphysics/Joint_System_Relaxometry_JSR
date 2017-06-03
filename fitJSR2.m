%% Function that estimates M0,T1,T2,B1 and B0 based on the JSR approach
% Usage:
%   [E_parameters, resNorm, res] = fitJSR(noisySignal,TRmprage,...
%     rageBlockDur,rageDeadTimeDur,rageTiDur,TRspgr,alpha_mprage,...
%     alpha_spgr,TRssfp,TE,alpha_ssfp,rf_phase_incr,rf_trms,x0,options)
%
%   noisySignal     -> measured signal - Expects magnitude MPRAGE and SPGR
%                       measurements as well as real and imaginary part of
%                       SSFP - [|MPRAGE(:)|; |SPGR(:)|; real(SSFP(:);
%                       imag(SSFP(:))];
%
%   rf_phase_incr   -> SSFP excitation phase increment at each RF in rad,
%                       can be an array of values
%
%   TRmprage        -> repetition time between MPRAGE 'alpha pulses' - RAGE
%                       block pulses time interval in ms, can be an array
%                       of values
%
%   rageBlockDur    -> Duration of the RAGE block = TFEfactor * TRmprage,
%                       can be an array of values
%
%   rageDeadTimeDur -> Duration between the last 'alpha pulse' of the RAGE
%                       block and the next inversion pulse in ms, can be an
%                       array of values
%
%   rageTiDur       -> Duration between the inversion pulse and the start
%                       of the RAGE block in ms
%
%   TRspgr          -> SPGR repetition time in ms, can be an array of
%                       values for each measurement
%
%   TRssfp          -> SSFP repetition time in ms, can be an array of
%                       values for each measurement
%
%   TE              -> readout echo time in ms, can be a single value which
%                       is applied to all measurements (recomended) but
%                       also accepts arreay of values with different
%                       readout durations for each separate measurement
%
%   alpha_mprage    -> RAGE block excitation pulses in deg, can accept an
%                       array of values corresponding to a different
%                       excitation per measrument
%
%   alpha_spgr      -> SPGR excitation pulses in deg, can accept an array
%                       of values corresponding to a different excitation
%                       per measurement
%
%   alpha_ssfp      -> SSFP excitation pulses in deg, can accept an array
%                       of values corresponding to a different excitation
%                       per measuremnt
%
%   rf_trms         -> Duration of SSFP RF pulse in ms
%
%   x0              -> Inital point of the fitting routine, can be empty
%                       and default value of [1 0 900 50 0 1] will be used
%
%   options         -> lsqnonlin object structure, can be empty and default
%                       structure will be used:
%                           options             = optimset('lsqnonlin');
%                           options.FinDiffType = 'central';
%                           options.TolFun      = 1e-12;
%                           options.TolX        = 1e-12;
%                           options.MaxFunEvals = 1e12;
%                           options.MaxIter     = 1e12;
%                           options.Display     = 'off';
%
%
% Rui Pedro A. G. Teixeira @ 31/05/2016
% Rui Pedro A. G. Teixeira @ 30/06/2016 - allow bounds to be set externaly
% Rui Pedro A. G. Teixeira @ 03/06/2017 - estimation of B1 making use of MPRAGE is being added but not functional yet it's WiP

function [E_parameters, resNorm, res, exitFlag] = fitJSR2(noisySignal,TRmprage,...
    rageBlockDur,rageDeadTimeDur,rageTiDur,TRspgr,alpha_mprage,...
    alpha_spgr,TRssfp,TE,alpha_ssfp,rf_phase_incr,rf_trms,x0,options,lb,ub,k)

%%
if ~isempty(k)
    fitb1 = false;
else
    fitb1 = true;
end

% Images are stored as single point precision but fit routines require
% double point precision
noisySignal     = double(noisySignal);

%% Scale the problem first by fitting M0
% This makes the problem easier to fit
% Also, check if we aim to estimate b1 or not and define cost function
% accordingly
if fitb1
    %% if user doesn't give a starting point we do an 'educated guess'
    if isempty(x0)
        x0 = [1 0 900 50 0 1];
        
    end
    SignalModel1     = @(E) expandToReal(-complex(E(1),E(2)),x0(3),x0(4),...
        -x0(5),rf_phase_incr,TRmprage,rageBlockDur,rageDeadTimeDur,...
        rageTiDur,TRspgr,TRssfp,TE,alpha_mprage,alpha_spgr,...
        alpha_ssfp,rf_trms,x0(6));
    
    CostFunction1    = @(E) noisySignal(:) - SignalModel1(E);
    
    [M0tmp,~,~,~] = lsqnonlin(CostFunction1,x0(1,1:2),[],[],options);
    
    % Now once a good starting point is defined do the actual full fit
    % forward model of the signal which is already expanded to be real and
    % imaginary separate as two different measurements -> avoids handling
    % complex data in lsqnonlin
    
    SignalModel     = @(E) expandToReal(-complex(E(1),E(2)),E(3),E(4),...
        -E(5),rf_phase_incr,TRmprage,rageBlockDur,rageDeadTimeDur,...
        rageTiDur,TRspgr,TRssfp,TE,alpha_mprage,alpha_spgr,...
        alpha_ssfp,rf_trms,E(6));
    
    % Cost Function is simply the difference between measurements and forward
    % model
    CostFunction    = @(E) ((noisySignal(:) - SignalModel(E)));
else
    %% if user doesn't give a starting point we do an 'educated guess'
    if isempty(x0)
        x0 = [1 0 900 50 0];
        
    end
%     
SignalModel1 = @(E) simulateDESPOTvoxel_complexCF(E(1)+1i.*E(2),x0(3),(x0(4)),(x0(5)),rf_phase_incr,...
            TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms);
        
    CostFunction1    = @(E) noisySignal(:) - SignalModel1(E);
%     
    [M0tmp,~,~,~] = lsqnonlin(CostFunction1,x0(1,1:2),[],[],options);
    
    % Now once a good starting point is defined do the actual full fit
    % forward model of the signal which is already expanded to be real and
    % imaginary separate as two different measurements -> avoids handling
    % complex data in lsqnonlin
    
    SignalModel = @(E) simulateJSRvoxel_complexCF(E(1)+1i.*E(2),(E(3)),(E(4)),(E(5)),rf_phase_incr,...
            TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms);
    
    % Cost Function is simply the difference between measurements and forward
    % model
    CostFunction    = @(E) ((noisySignal(:) - SignalModel(E)));
end


% if options are not given as input create option structure
if isempty(options)
    
    options             = optimset('lsqnonlin');
    options.FinDiffType = 'central';
    options.TolFun      = 1e-6;
    options.TolX        = 1e-6;
    options.MaxFunEvals = 1e12;
    options.MaxIter     = 1e2;
    options.Display     = 'off';
    
end
%%
x0 = [M0tmp(:);x0(3:end)'];
%%
[E_parameters,resNorm,res,exitFlag] = lsqnonlin(CostFunction,x0,lb,ub,options);
end
