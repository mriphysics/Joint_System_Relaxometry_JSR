%%
clear;clc;close all;
M_0i = 0;
M_0r = 10;
RF_dur = 0.6;
T_1 = 1000;
T_2 = 60;
T_E = 2;
T_R_spgr = 8.2;
T_R_ssfp = 4;
a_spgr = deg2rad([6 8 10 12 14 16]);
a_ssfp = deg2rad([15 25 35 45 55 65 25 55]);
b_0 = 10*2*pi*T_R_ssfp*10^-3;
phi_incr = [pi pi pi pi pi pi 0 0];
tic
% for ii=1:1e2
tmp = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr(:) );
% end
toc
figure(1);
plot([a_spgr],tmp(1:length(a_spgr)));
hold all;
plot([a_ssfp],tmp(length(a_spgr)+1:length(a_spgr)+length(a_ssfp)));
plot([a_ssfp],tmp(length(a_spgr)+length(a_ssfp)+1:length(a_spgr)+length(a_ssfp)+length(a_ssfp)));
hold off;
xlabel('Flip angle (rad)')
ylabel('Signal Intensity')
%%
[tmp,tmpJacobian] = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr(:) );
F = ((tmpJacobian.'/(diag(0.02*ones(length([a_spgr a_ssfp a_ssfp]),1),0))*tmpJacobian));
CRLB = sqrt(diag(inv(F)))
%%
clear

M_0i = 0;
M_0r = 10;
T_1 = 2063;
T_2 = 184;
RF_dur = 0.6;
% [T1, T2] = ndgrid(600:25:1200, 25:5:80);

sz = size(T_1);
k = 1;
% b_0 = 10*2*pi*T_R_ssfp;

T_R_spgr = 6.2;

a_spgr = deg2rad([4 18]);

TRmprage = [];
alpha_mprage = [];
rageBlockDur = [];
rageTiDur = [];
rageDeadTimeDur = [] ;

T_R_ssfp = 4.2;
b_0 = 10*2*pi*T_R_ssfp*10^-3;
a_ssfp = deg2rad([15 65 15 65]);
phi_incr = [pi pi 0 0];

T_E = 2.1;
RF_dur = 0.6;
Ntrials = 1e3;
stdMatrix = zeros(sz(1),sz(2),5);
crlbMatrix = zeros(sz(1),sz(2),5);
%
options                 = optimset('lsqnonlin');
options.FinDiffType     = 'forward';
options.Algorithm = {'levenberg-marquardt',1e6};
% options.PlotFcns = {@optimplotresnorm,@optimplotx};
options.TolFun          = 1e-15;
options.TolX            = 1e-15;
options.MaxFunEvals     = 1e12;
options.MaxIter         = 1e12;
options.Display         = 'off';
options.Jacobian        = 'on';
%
% profile on
tic
for T1idx = 1:sz(1)
    for T2idx = 1:sz(2)
        
        %%
        [model] = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr );
        
        stdNoise = 0.002*M_0r;
        noisyModel = repmat(model(:)',Ntrials,1) + stdNoise.*randn(Ntrials,10);
%         noisyModel(:,1:2) = abs(noisyModel(:,1:2));
        E_parameters = zeros(Ntrials,5);
%         E_DESPOT1    = zeros(Ntrials,2);
%         flags.k = k;
%         flags.fitb1 = false;
        
        for ii=1:Ntrials
            [E_parameters(ii,:)] = fit_jsr_model(noisyModel(ii,:),...
        RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr,...
        [],options,[],[]);
%             [E_DESPOT1(ii,:)] = fitDESPOT1(noisyModel(ii,1:2),TRspgr,alpha_spgr);
            %     disp(E_parameters(ii,:))
        end
%         stdMatrix(T1idx,T2idx,:) = std(E_parameters);
%         crlbMatrix(T1idx,T2idx,:) = CRLB_JSR_noB1(M0,T1(T1idx,T2idx),T2(T1idx,T2idx),b0,rf_phase_incr,...
%             TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms,k);
        disp(['T1idx = ' num2str(T1idx) ' T2idx = ' num2str(T2idx)])
    end
%     disp(T1idx)
%     toc
%     save('MonteCarloCRLB_Validation_noB1');
end
toc
disp(E_parameters(1,:));
% profile viewer
% profile off
% save('MonteCarloCRLB_Validation_noB1');
%
figure(1);
histogram(E_parameters(:,3),500:50:3500);
% hold all;
% histogram(E_DESPOT1(:,2),500:50:3500);
% hold off;