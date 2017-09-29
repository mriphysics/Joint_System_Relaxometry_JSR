%% Find optimal FA set
function [optimizedFAset,fval] = JSRplus_FindOptimSet(rM0grid, iM0grid, R1grid, R2grid, b0grid, b1grid,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,varargin)
%%

minTR   = 7;%min TRssfp = 7ms;
maxTR   = 100;
Ttotal  = 96;

% Number of mprage measurements
Nmprage = length(alpha_mprage);

% Number of mprage measurements
Nspgr = length(alpha_spgr);

% Number of ssfp measruments
Nssfp = length(alpha_ssfp);


%% Parameters to optimize:
% I'll start by defining the indexes of each parameter - makes code more
% readable.

%
if (Nmprage > 0)&&(Nspgr > 0)
    %% Start with the case where all sequences are used:
    TRmprageIdx     = 1                         : Nmprage;
    alpha_mprageIdx = TRmprageIdx(end)      + 1 : TRmprageIdx(end)      + Nmprage;
    TFE_factorIdx   = alpha_mprageIdx(end)  + 1 : alpha_mprageIdx(end)  + Nmprage;
    
    TRspgrIdx       = TFE_factorIdx(end)    + 1 : TFE_factorIdx(end)    + Nspgr;
    alpha_spgrIdx   = TRspgrIdx(end)        + 1 : TRspgrIdx(end)        + Nspgr;
    
    alpha_ssfpIdx   = alpha_spgrIdx(end)    + 1 : alpha_spgrIdx(end)    + Nssfp;
    phi_ssfpIdx     = alpha_ssfpIdx(end)    + 1 : alpha_ssfpIdx(end)    + Nssfp;
    
    % the non linear constraint needs to be dynamically generated
    disp('Generating constraint file...');
    
    FID = fopen('totalTimeConstraint.m','w+');
    fprintf(FID,'function [c, ceq] = totalTimeConstraint(x)\n');
    
    fprintf(FID,['c = [sum(x([' repmat('%d ',1,Nmprage + Nspgr) '])) + %d .* %f - %f];\n'],[[TRmprageIdx TRspgrIdx]  Nssfp minTR Ttotal]);
    fprintf(FID,'ceq = [];\n');
    fprintf(FID,'end');
    fclose(FID);
    disp('Generating constraint file... Done');
    
    nonlcon = @totalTimeConstraint;
    
    % Define the starting point for the optimization
    x0(TRmprageIdx)     = TRmprage;
    x0(alpha_mprageIdx) = alpha_mprage;
    x0(TFE_factorIdx)   = t_rage./TRssfp;
    x0(TRspgrIdx)       = TRspgr;
    x0(alpha_spgrIdx)   = alpha_spgr;
    x0(alpha_ssfpIdx)   = alpha_ssfp;
    x0(phi_ssfpIdx)     = phi_ssfp;
    
elseif (Nmprage == 0)&&(Nspgr > 0)
    
    %% Next step consider when there's no MPRAGEand all other sequences are used:
    TRmprageIdx     = [];
    alpha_mprageIdx = [];
    TFE_factorIdx   = [];
    
    TRspgrIdx       = 1 : Nspgr;
    alpha_spgrIdx   = TRspgrIdx(end)        + 1 : TRspgrIdx(end)        + Nspgr;
    
    alpha_ssfpIdx   = alpha_spgrIdx(end)    + 1 : alpha_spgrIdx(end)    + Nssfp;
    phi_ssfpIdx     = alpha_ssfpIdx(end)    + 1 : alpha_ssfpIdx(end)    + Nssfp;
    
    % the non linear constraint needs to be dynamically generated
    disp('Generating constraint file...');
    
    FID = fopen('totalTimeConstraint.m','w+');
    fprintf(FID,'function [c, ceq] = totalTimeConstraint(x)\n');
    
    fprintf(FID,['c = [sum(x([' repmat('%d ',1,Nmprage + Nspgr) '])) + %d .* %f - %f];\n'],[[TRmprageIdx TRspgrIdx]  Nssfp minTR Ttotal]);
    fprintf(FID,'ceq = [];\n');
    fprintf(FID,'end');
    fclose(FID);
    disp('Generating constraint file... Done');
    
    nonlcon = @totalTimeConstraint;
    
    % Define the starting point for the optimization
    x0(TRmprageIdx)     = TRmprage;
    x0(alpha_mprageIdx) = alpha_mprage;
    x0(TFE_factorIdx)   = t_rage./TRssfp;
    x0(TRspgrIdx)       = TRspgr;
    x0(alpha_spgrIdx)   = alpha_spgr;
    x0(alpha_ssfpIdx)   = alpha_ssfp;
    x0(phi_ssfpIdx)     = phi_ssfp;
    
elseif (Nmprage > 0)&&(Nspgr == 0)
    
    %% Next step consider when there's no SPGR and all other sequences are used:
    TRmprageIdx     = 1                         : Nmprage;
    alpha_mprageIdx = TRmprageIdx(end)      + 1 : TRmprageIdx(end)      + Nmprage;
    TFE_factorIdx   = alpha_mprageIdx(end)  + 1 : alpha_mprageIdx(end)  + Nmprage;
    
    TRspgrIdx       = [];
    alpha_spgrIdx   = [];
    
    alpha_ssfpIdx   = TRmprageIdx(end)      + 1 : TRmprageIdx(end)    + Nssfp;
    phi_ssfpIdx     = alpha_ssfpIdx(end)    + 1 : alpha_ssfpIdx(end)    + Nssfp;
    
    % the non linear constraint needs to be dynamically generated
    disp('Generating constraint file...');
    
    FID = fopen('totalTimeConstraint.m','w+');
    fprintf(FID,'function [c, ceq] = totalTimeConstraint(x)\n');
    
    fprintf(FID,['c = [sum(x([' repmat('%d ',1,Nmprage + Nspgr) '])) + %d .* %f - %f];\n'],[[TRmprageIdx TRspgrIdx]  Nssfp minTR Ttotal]);
    fprintf(FID,'ceq = [];\n');
    fprintf(FID,'end');
    fclose(FID);
    disp('Generating constraint file... Done');
    
    nonlcon = @totalTimeConstraint;
    
    % Define the starting point for the optimization
    x0(TRmprageIdx)     = TRmprage;
    x0(alpha_mprageIdx) = alpha_mprage;
    x0(TFE_factorIdx)   = t_rage./TRssfp;
    x0(TRspgrIdx)       = TRspgr;
    x0(alpha_spgrIdx)   = alpha_spgr;
    x0(alpha_ssfpIdx)   = alpha_ssfp;
    x0(phi_ssfpIdx)     = phi_ssfp;
    
elseif (Nmprage == 0)&&(Nspgr == 0)
    %% Next step consider when there's no SPGR and all other sequences are used:
    TRmprageIdx     = [];
    alpha_mprageIdx = [];
    TFE_factorIdx   = [];
    
    TRspgrIdx       = [];
    alpha_spgrIdx   = [];
    
    alpha_ssfpIdx   = 1                         : Nssfp;
    phi_ssfpIdx     = alpha_ssfpIdx(end)    + 1 : alpha_ssfpIdx(end)  + Nssfp;
    
    % the non linear constraint needs to be dynamically generated
    disp('Generating constraint file...');
    
    FID = fopen('totalTimeConstraint.m','w+');
    fprintf(FID,'function [c, ceq] = totalTimeConstraint(x)\n');
    
    fprintf(FID,['c = [sum(x([' repmat('%d ',1,Nmprage + Nspgr) '])) + %d .* %f - %f];\n'],[[TRmprageIdx TRspgrIdx]  Nssfp minTR Ttotal]);
    fprintf(FID,'ceq = [];\n');
    fprintf(FID,'end');
    fclose(FID);
    disp('Generating constraint file... Done');
    
    nonlcon = @totalTimeConstraint;
    
    % Define the starting point for the optimization
    x0(TRmprageIdx)     = TRmprage;
    x0(alpha_mprageIdx) = alpha_mprage;
    x0(TFE_factorIdx)   = t_rage./TRssfp;
    x0(TRspgrIdx)       = TRspgr;
    x0(alpha_spgrIdx)   = alpha_spgr;
    x0(alpha_ssfpIdx)   = alpha_ssfp;
    x0(phi_ssfpIdx)     = phi_ssfp;
    
end


%% Now that the indexes are created it becomes simple to create the bounds
% Define lower bounds for parameters to be optimized

lb(TRmprageIdx)     = 12; % MPRAGE minTR is different due to SAR constraints of 180 pulse
lb(alpha_mprageIdx) = deg2rad(6);
lb(TFE_factorIdx)   = 50;
lb(TRspgrIdx)       = minTR;
lb(alpha_spgrIdx)   = deg2rad(6);
lb(alpha_ssfpIdx)   = deg2rad(6);
lb(phi_ssfpIdx)     = phi_ssfp;

% Define upper bounds for parameters to be optimized
ub(TRmprageIdx)     = maxTR;
ub(alpha_mprageIdx) = deg2rad(90);
ub(TFE_factorIdx)   = 1000;
ub(TRspgrIdx)       = maxTR;
ub(alpha_spgrIdx)   = deg2rad(90);
ub(alpha_ssfpIdx)   = deg2rad(90);
ub(phi_ssfpIdx)     = phi_ssfp;


%% Create cost function

CostFunction = @(X) JSRplus_CRLB_grid(rM0grid, iM0grid, R1grid, R2grid, b0grid, b1grid,...
    X(TRmprageIdx), TEmprage, X(alpha_mprageIdx), X(TRmprageIdx).*round(X(TFE_factorIdx)), X(TRmprageIdx), X(TRmprageIdx), ...
    X(TRspgrIdx), TEspgr, X(alpha_spgrIdx),...
    TRssfp, TEssfp, X(alpha_ssfpIdx), X(phi_ssfpIdx), rf_trms);


%%
A   = [];
b   = [];
Aeq = [];
beq = [];
% lb = [];
% ub = [];
% nonlcon = @mprageConstraint;
options = psoptimset('Display','iter','PlotFcns',{@psplotbestf, @psplotbestx, @psplotmeshsize});
options.CompletePoll = 'on';
options.CompleteSearch = 'on';
options.MeshAccelerator = 'on';
options.Cache = 'on';
options.TolFun = 1e-1;
options.TolMesh = 1e-3;
% options.InitialMeshSize = 1e3;
options.MaxIter = 5;
options.MaxFunEvals = 1e12;
% options.InitialPenalty = 1e4;
options.ScaleMesh = 'off';
% options.PollMethod = 'random';
options.UseParallel = false;
% options.MaxIter     = 310;

try
    [optimizedFAset,fval] = patternsearch(@(X) 100*CostFunction(X),x0(:),A,b,Aeq,beq,lb(:),ub(:),nonlcon,options);
catch ME
    warning('Problem using function.  Assigning a value of []');
    optimizedFAset = [];
    fval    = [];
end
