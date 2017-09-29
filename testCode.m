%% Script which seeks to test code for JSR+
% tidy up workspace and command window
clear;
clc;

% generate relaxometry parameters
M0 = complex(10,10);
R1 = 1/1000;
R2 = 1/50;
b0 = 0;
b1 = 1;

% Create Acquistion parameters
TRmprage = 12;
TEmprage = 6.0;
t_rage   = 1212;
TI = 12;
TD = 12;
alpha_mprage = deg2rad([6.000,  10.000,  14.000,  18.000,  22.000])';

% TRspgr = 12;
% TEspgr = 6;
% alpha_spgr = deg2rad([6.000,  10.000,  14.000,  18.000,  22.000])';

TRspgr = [];
TEspgr = [];
alpha_spgr = [];


TRssfp = 12;
TEssfp = 6;
alpha_ssfp = (deg2rad([15.000,  30.000,  45.000,  60.000,  75.000,  15.000,  30.000,  45.000,  60.000,  75.000])');
phi_ssfp = deg2rad([180, 180, 180, 180, 180,   0,   0,   0,   0,   0]');
rf_trms = (1.2*3.0/2.1362);

% Simulate actual sequences
[S,J] = simulateJSRplus(M0, R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);

% plor data points
figure(1);
plot(rad2deg([alpha_mprage;alpha_spgr;alpha_ssfp;alpha_ssfp]),S,'o')

% Compute CRLB of JSR+
noise_var = 0.002;
covM = noise_var*eye(length([alpha_mprage;alpha_spgr;alpha_ssfp;alpha_ssfp]));
F = J'/covM*J;
CRLB = sqrt(diag(F\eye(size(F))));
%%  Now proceed with a MonteCarlo trial
% here I proceed not only the ability to estimate all the paramters but
% also if the CRLB allows us to estimate the variance of the estimation in
% other words, the good old sanity check
% Start, for simplicity with everything having gaussian distributed noise
Ntrials = 1e3;

noisyS = repmat(S,1,Ntrials) + sqrt(noise_var).*randn(length(S),Ntrials);
noisyS = noisyS';

E_par = zeros(Ntrials,6);
outputFlag = zeros(Ntrials,1);
tic
options             = optimset('lsqnonlin');
options.Algorithm   = 'levenberg-marquardt';
options.TolFun      = 1e-12;
options.TolX        = 1e-12;
options.MaxIter     = 50;
options.MaxFunEvals = 1e3;
options.Display     = 'off';
options.Jacobian    = 'on';
options.CheckGradients = false;
options.ScaleProblem = 'jacobian';
parfor ii=1:Ntrials
    [E_par(ii,:), ~, ~, outputFlag(ii)] = fitJSRplus(noisyS(ii,:),...
        TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
        TRspgr, TEspgr, alpha_spgr,...
        TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,...
        [real(M0),imag(M0),R1,R2,b0,b1],options,[],[]);
end
toc
%%
figure(2);
GroundTruth = [real(M0) imag(M0) R1 R2 b0 b1];
titleLabels = {'rM_0','iM_0','R_1','R_2','B_0','B_1'};
std(E_par);
for ii=1:6
    subplot(2,3,ii);
    histogram(E_par(:,ii),'Normalization','pdf');
    hold all;
    yl = ylim;
    xl = xlim;
    plot([GroundTruth(ii) GroundTruth(ii)],[yl(1) yl(2)],'LineWidth',2);
    hold off;
    title(titleLabels{ii},'FontSize',20,'FontWeight','bold')
end

MC_STD = std(E_par)';
VarName = {'rM0','iM0','R_1','R_2','B_0','B_1'};
T = table(CRLB,MC_STD,'RowNames',VarName);
disp(T);
% [CRLB std(E_par)']
%%
