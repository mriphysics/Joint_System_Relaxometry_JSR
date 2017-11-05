%% Script which seeks to test code for JSR+
% tidy up workspace and command window
% clear;
% clc;

% generate relaxometry parameters
M0 = complex(-900,-300)./100;
R1 = 1200;
R2 = 30;
b0 = 0;
b1 = 1;

% Create Acquistion parameters
TRmprage = 14;
TEmprage = 3.5;
t_rage   = [236*TRmprage]';
TI = TRmprage;
TD = TRmprage;
alpha_mprage = deg2rad([6.8])';

% TRmprage = [];
% TEmprage = [];
% t_rage   = [];
% TI = [];
% TD = [];
% alpha_mprage = deg2rad([])';

TRspgr = [14]';
TEspgr = 3.5;
alpha_spgr = deg2rad([15.9])';

% TRspgr = [];
% TEspgr = [];
% alpha_spgr = [];


TRssfp = 7;
TEssfp = 3.5;
% Adult WM optimized Protocol
alpha_ssfp = deg2rad([31.7 6.5 20.6 14])'; 
phi_ssfp = deg2rad([290 298 39 223]');
% Phantom / Neonatal WM optimized Protocol
% alpha_ssfp = deg2rad([10 12 17.7 9.5 90 21.2 13.3 90 10])'; 
% phi_ssfp = deg2rad(round([189 351 345 6 180 18 43 180 189]'));
rf_trms = (1.2*3.0/2.1362);

% Simulate actual sequences
[S,J] = simulateJSRplus(M0, R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
% [S] = simulateJSRplus(M0, R1, R2, b0, b1,...
%     TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
%     TRspgr, TEspgr, alpha_spgr,...
%     TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
% plor data points
figure(42);
plot(rad2deg([alpha_mprage;alpha_spgr;alpha_ssfp;alpha_ssfp]),S,'o')

% Compute CRLB of JSR+
noise_var = 0.0008;
SNR = 0.1*abs(M0)/sqrt(noise_var);
disp(['SNR: ' num2str(SNR)]);
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
options.Jacobian    = 'off';
options.CheckGradients = false;
options.ScaleProblem = 'jacobian';
options(2) = options;
parfor ii=1:Ntrials
    [E_par(ii,:), ~, ~, outputFlag(ii)] = fitJSRplus(noisyS(ii,:),...
        TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
        TRspgr, TEspgr, alpha_spgr,...
        TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,...
        [10 0 1000 50 0 1],options,[],[]);
end
toc
%%
figure(2);
GroundTruth = [real(M0) imag(M0) R1 R2 b0 b1]';
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
T = table(CRLB,MC_STD,GroundTruth,'RowNames',VarName);
disp(T);
% [CRLB std(E_par)']

%%
% tic
% [R1grid, R2grid] = ndgrid(0:1e-5:4e-2,0:1e-5:5e-2);
% CF_MPRAGE = zeros(length(R1grid(:)),1);
% CF_noMPRAGE = CF_MPRAGE;
% for ii=1:length(R1grid(:))
%     
% [Snew_MPRAGE] = simulateJSRplus(M0, R1grid(ii), R2grid(ii), b0, b1,...
%     TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
%     TRspgr, TEspgr, alpha_spgr,...
%     TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
% [Snew_noMPRAGE] = simulateJSRplus(M0, R1grid(ii), R2grid(ii), b0, b1,...
%     [], [], [], [], TI, TD, ...
%     TRspgr, TEspgr, alpha_spgr,...
%     TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
% 
% CF_MPRAGE(ii) = sqrt(sum(S(:).^2 + Snew_MPRAGE(:).^2));
% CF_noMPRAGE(ii) = sqrt(sum(S(2:end).^2 + Snew_noMPRAGE(:).^2));
% end
% CF_MPRAGE = reshape(CF_MPRAGE,size(R1grid,1),size(R1grid,2));
% CF_noMPRAGE = reshape(CF_noMPRAGE,size(R1grid,1),size(R1grid,2));
% toc
% %%
% figure;
% subplot(1,2,1);
% imagesc(R2grid(1,:),R1grid(:,1),CF_MPRAGE,[0 40]);
% 
% subplot(1,2,2);
% imagesc(CF_noMPRAGE,[0 40]);