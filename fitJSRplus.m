%% Code Generated to estimate JSR+ parameters

function [E_parameters, resNorm, res, exitFlag] = fitJSRplus(noisyS,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,...
    x0,options,lb,ub,varargin)

% Let's try and make the code general I'll set some default flags that can
% be changed depending on which options we wish to build

knownB1 = false; % default is that we don't know B1 and which to estimate it

for ii=1:length(varargin)
    if strcmp(varargin{ii},'knownB1')
        knownB1 = true;
        b1 = varargin{ii+1};
    end
end

% If starting point is not provided make educated guess of white matter
if isempty(x0)
    x0 = [10 0 1/700 1/60 0 1];
end

% If fitting options are not provided go with default settings
if isempty(options)
    options             = optimset('lsqnonlin');
    options.Algorithm   = 'levenberg-marquardt';
    options.TolFun      = 1e-6;
    options.TolX        = 1e-6;
    options.MaxIter     = 50;
    options.Display     = 'off';
    
    % LM algorithm does not accept bounds
    lb                  = [];
    ub                  = [];
end

if knownB1
    try
        [S] = @(E) simulateJSRplus(E(1)+1i*E(2), x0(3), x0(4), x0(5), b1,...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        CF = @(E) (S(E) - noisyS(:));
        [E_tmp] = lsqnonlin(CF,x0(1:2),lb,ub,options(1));
        x0(1:2) = E_tmp;
        x0(6) = [];
        
        [S] = @(E) simulateJSRplus(E(1)+1i*E(2), E(3), E(4), E(5), b1,...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        CF = @(E) (S(E) - noisyS(:));
        [E_parameters,resNorm,res,exitFlag] = lsqnonlin(CF,x0,lb,ub,options(2));
    catch ME
        E_parameters = zeros(1,5);
        resNorm = 0;
        res = zeros(1,length([alpha_mprage(:); alpha_spgr(:); alpha_ssfp(:);alpha_ssfp(:)]));
        exitFlag = 999;
        disp('Error in lsqnonlin');
    end
else
    try
        [S] = @(E) simulateJSRplus(E(1)+1i*E(2), x0(3), x0(4), x0(5), x0(6),...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        CF = @(E) (S(E) - noisyS(:));
        [E_tmp] = lsqnonlin(CF,x0(1:2),lb,ub,options(1));
        x0(1:2) = E_tmp;
        
        [S] = @(E) simulateJSRplus(E(1)+1i*E(2), E(3), E(4), E(5), E(6),...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        CF = @(E) (S(E) - noisyS(:));
        [E_parameters,resNorm,res,exitFlag] = lsqnonlin(CF,x0,lb,ub,options(2));
    catch ME
        E_parameters = zeros(1,6);
        resNorm = 0;
        res = zeros(1,length([alpha_mprage(:); alpha_spgr(:); alpha_ssfp(:);alpha_ssfp(:)]));
        exitFlag = 999;
        disp('Error in lsqnonlin');
    end
end

function [CF,J_CF] = JSRplus_CF(E,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,noisyS)

[S,J_CF] = simulateJSRplus(E(1)+1i*E(2), E(3), E(4), E(5), E(6),...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);

CF = S(:) - noisyS(:);