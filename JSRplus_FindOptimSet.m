%% Find optimal FA set
function [optimizedFAset,fval] = JSRplus_FindOptimSet(rM0grid, iM0grid, R1grid, R2grid, b0grid, b1grid,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,varargin)
%%
knownB1 = false;
for ii=1:length(varargin)
    if strcmp(varargin{ii},'knownB1')
        knownB1 = true;
    end
end

minTR   = 7;%min TRssfp = 7ms;
maxTR   = 20;
Ttotal  = 57;

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
    x0(TFE_factorIdx)   = round(t_rage./TRmprage)/1000;
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
    x0(TFE_factorIdx)   = round(t_rage./TRmprage)/1000;
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
    
    alpha_ssfpIdx   = TFE_factorIdx(end)      + 1 : TFE_factorIdx(end)    + Nssfp;
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
    x0(TFE_factorIdx)   = round(t_rage./TRmprage)/1000;
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
    x0(TFE_factorIdx)   = round(t_rage./TRmprage)/1000;
    x0(TRspgrIdx)       = TRspgr;
    x0(alpha_spgrIdx)   = alpha_spgr;
    x0(alpha_ssfpIdx)   = alpha_ssfp;
    x0(phi_ssfpIdx)     = phi_ssfp;
    
end


%% Now that the indexes are created it becomes simple to create the bounds
% Define lower bounds for parameters to be optimized

lb(TRmprageIdx)     = 14;%12; % MPRAGE minTR is different due to SAR constraints of 180 pulse
lb(alpha_mprageIdx) = deg2rad(6);
lb(TFE_factorIdx)   = 200/100;
lb(TRspgrIdx)       = 14;
lb(alpha_spgrIdx)   = deg2rad(6);
lb(alpha_ssfpIdx)   = deg2rad(6);
lb(phi_ssfpIdx)     = 0;

% Define upper bounds for parameters to be optimized
ub(TRmprageIdx)     = 14;%15;
ub(alpha_mprageIdx) = deg2rad(90);
ub(TFE_factorIdx)   = 900/100;
ub(TRspgrIdx)       = 14;
ub(alpha_spgrIdx)   = deg2rad(90);
ub(alpha_ssfpIdx)   = deg2rad(90);
ub(phi_ssfpIdx)     = 2*pi;


%% Create cost function

CostFunction = @(X) JSRplus_CRLB_grid(rM0grid, iM0grid, R1grid, R2grid, b0grid, b1grid,...
    X(TRmprageIdx), TEmprage, X(alpha_mprageIdx), X(TRmprageIdx).*100.*round(X(TFE_factorIdx)), X(TRmprageIdx), X(TRmprageIdx), ...
    X(TRspgrIdx), TEspgr, X(alpha_spgrIdx),...
    TRssfp, TEssfp, X(alpha_ssfpIdx), X(phi_ssfpIdx), rf_trms,knownB1);

%%
UsePatternSearch = false;
UseSimAnnealing = true;
UseFmincon = false;
UseDirectSearch = false;
if UsePatternSearch
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    % lb = [];
    % ub = [];
    % nonlcon = @mprageConstraint;
    options = psoptimset('Display','final','PlotFcns',{@psplotbestf, @psplotbestx, @psplotmeshsize});
%     options.CompletePoll = 'on';
%     options.CompleteSearch = 'on';
%     options.MeshAccelerator = 'off';
%     options.Cache = 'on';
    options.TolFun = 1;
%     options.TolMesh = 1e-6;
    % options.InitialMeshSize = 1e3;
    options.MaxIter = 5;
    % options.MaxFunEvals = 1e12;
    % options.InitialPenalty = 1e4;
%     options.ScaleMesh = 'on';
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
elseif UseSimAnnealing
    
    options = saoptimset('PlotFcns',{@saplotbestx,...
        @saplotbestf,@saplotx,@saplotf});
    options.Display = 'final';
%     options.TolFun = 0.1;
%     options.HybridFcn = @fmincon;
%     options.HybridInterval = 'end';
%     options.ReannealInterval = 200;
%     options.InitialTemperature = 1e4;
    options.AnnealingFcn = @annealingboltz;
%     options.MaxIter = 2000;
%     options.MaxStallIterations = 3000;
    try
        [optimizedFAset,fval] = simulannealbnd(@(X) real(100*CostFunction(X)),x0(:),lb(:),ub(:),options);
    catch ME
        warning('Problem using function.  Assigning a value of []');
        optimizedFAset = [];
        fval    = [];
    end
end    
if UseFmincon
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    % lb = [];
    % ub = [];
    % nonlcon = @mprageConstraint;
    options = optimoptions('fmincon','Display','iter');%,'PlotFcns',{@optimplotx, @optimplotfunccount , @optimplotfval ,@optimplotconstrviolation ,optimplotstepsize ,optimplotfirstorderopt });
    options.PlotFcn = {@optimplotx, @optimplotfunccount , @optimplotfval ,@optimplotconstrviolation ,@optimplotstepsize ,@optimplotfirstorderopt };
    
    try
        [optimizedFAset,fval] = fmincon(@(X) 100*CostFunction(X),optimizedFAset(:),A,b,Aeq,beq,lb(:),ub(:),nonlcon,options);
    catch ME
        warning('Problem using function.  Assigning a value of []');
        optimizedFAset = [];
        fval    = [];
    end
elseif UseDirectSearch
    CF = @(X) real(100*CostFunction(X));
    stepGrid = 0.1;
    TRmprageGrid = lb(TRmprageIdx(1)):stepGrid*(ub(TRmprageIdx(1)) - lb(TRmprageIdx(1))):ub(TRmprageIdx(1));
    TRmprageGrid = repmat(TRmprageGrid(:),1,length(TRmprageIdx));

    alpha_mprageGrid = lb(alpha_mprageIdx(1)):stepGrid*(ub(alpha_mprageIdx(1)) - lb(alpha_mprageIdx(1))):ub(alpha_mprageIdx(1));
    alpha_mprageGrid = repmat(alpha_mprageGrid(:),1,length(alpha_mprageIdx));

    TFE_factorGrid = lb(TFE_factorIdx(1)):stepGrid*(ub(TFE_factorIdx(1)) - lb(TFE_factorIdx(1))):ub(TFE_factorIdx(1));
    TFE_factorGrid = repmat(TFE_factorGrid(:),1,length(TFE_factorIdx));
    
    TRspgrGrid = lb(TRspgrIdx(1)):stepGrid*(ub(TRspgrIdx(1)) - lb(TRspgrIdx(1))):ub(TRspgrIdx(1));
    TRspgrGrid = repmat(TRspgrGrid(:),1,length(TRspgrIdx));

    alpha_spgrGrid = lb(alpha_spgrIdx(1)):stepGrid*(ub(alpha_spgrIdx(1)) - lb(alpha_spgrIdx(1))):ub(alpha_spgrIdx(1));
    alpha_spgrGrid = repmat(alpha_spgrGrid(:),1,length(alpha_spgrIdx));
    
    alpha_ssfpGrid = lb(alpha_ssfpIdx(1)):stepGrid*(ub(alpha_ssfpIdx(1)) - lb(alpha_ssfpIdx(1))):ub(alpha_ssfpIdx(1));
    alpha_ssfpGrid = repmat(alpha_ssfpGrid(:),1,length(alpha_ssfpIdx));

    phi_ssfpGrid = lb(phi_ssfpIdx(1)):stepGrid*(ub(phi_ssfpIdx(1)) - lb(phi_ssfpIdx(1))):ub(phi_ssfpIdx(1));
    phi_ssfpGrid = repmat(phi_ssfpGrid(:),1,length(phi_ssfpIdx));
    optimizedFAset = [];
    fval = [];
    BestCF = 1e4;
    
for aa = 1:size(TRmprageGrid,1)
    for bb = 1:size(alpha_mprageGrid,1)
        for cc = 1:size(TFE_factorGrid,1)
            for dd = 1:size(TRspgrGrid,1)
                for ee = 1:size(alpha_spgrGrid,1)
                    for ff = 1:size(alpha_ssfpGrid,1)
                        for gg = 1:size(phi_ssfpGrid,1)
                            protocol = [TRmprageGrid(aa,:) alpha_mprageGrid(bb,:) TFE_factorGrid(cc,:) TRspgrGrid(dd,:) alpha_spgrGrid(ee,:) alpha_ssfpGrid(ff,:) phi_ssfpGrid(gg,:)];
                            tmp = CF(protocol);
                            if tmp<BestCF
                                fval = tmp;
                                optimizedFAset = protocol;
                            end
                        end
                    end
                end
            end
        end
    end
end
end


