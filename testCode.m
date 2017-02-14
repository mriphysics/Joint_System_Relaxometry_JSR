%% Test script of JSR code
% this script allows a quick example of how to use the functions in this code

%% As a quick example here it is shown how to generate the signal curves

clear;clc;close all;
% M0 = M_0r + 1iM_0i
M_0i = 0;
M_0r = 10;
% rf pulse duration
RF_dur = 0.6;
T_1 = 1000;
T_2 = 60;
T_E = 2;
T_R_spgr = 8.2;
T_R_ssfp = 4;
a_spgr = deg2rad([6 8 10 12 14 16]);
a_ssfp = deg2rad([15 25 35 45 55 65 25 55]);
b_0 = 10*2*pi*T_R_ssfp*10^-3; % b_0 is defined has dephasing per ssfp repetition time, in this example 10 is the underlying field map in Hz
phi_incr = [pi pi pi pi pi pi 0 0];

% This generates the actual jsr signal as a single array
tic
tmp = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr(:) );
toc

% plot the generate signal as a function of flip angle
figure(1);
plot([a_spgr],tmp(1:length(a_spgr)));
hold all;
plot([a_ssfp],tmp(length(a_spgr)+1:length(a_spgr)+length(a_ssfp)));
plot([a_ssfp],tmp(length(a_spgr)+length(a_ssfp)+1:length(a_spgr)+length(a_ssfp)+length(a_ssfp)));
hold off;
xlabel('Flip angle (rad)')
ylabel('Signal Intensity')

%% Example to export the jacobian of the signal at the given paramters
% The jacobian is important not only to make the fit more efficient it can be used to compute the CRLB
[tmp,tmpJacobian] = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr(:) );

F = ((tmpJacobian.'/(diag(0.02*ones(length([a_spgr a_ssfp a_ssfp]),1),0))*tmpJacobian)); % computes the Fisher information matrix
CRLB = sqrt(diag(inv(F))); % CRLB is the inverse of the Fisher matrix

%% Here we do a small monte carlo simulation to demonstrate how to apply the fitting routines
clear

M_0i = 0;
M_0r = 10;
T_1 = 2063;
T_2 = 184;
RF_dur = 0.6;

sz = size(T_1);
k = 1;

T_R_spgr = 6.2;

a_spgr = deg2rad([4 18]);

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
options.TolFun          = 1e-15;
options.TolX            = 1e-15;
options.MaxFunEvals     = 1e12;
options.MaxIter         = 1e12;
options.Display         = 'off';
options.Jacobian        = 'on';

% profile on
tic
for T1idx = 1:sz(1)
    for T2idx = 1:sz(2)
        
        %%
        [model] = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr );
        
        stdNoise = 0.002*M_0r;
        noisyModel = repmat(model(:)',Ntrials,1) + stdNoise.*randn(Ntrials,10);
        E_parameters = zeros(Ntrials,5);
        
        for ii=1:Ntrials
            [E_parameters(ii,:)] = fit_jsr_model(noisyModel(ii,:),...
        RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr,...
        [],options,[],[]);
        end
        disp(['T1idx = ' num2str(T1idx) ' T2idx = ' num2str(T2idx)])
    end
end
toc
disp(E_parameters(1,:));
% profile viewer
% profile off

% plot a histogram of the final Fit
figure(1);
subplot(1,2,1);
histogram(E_parameters(:,3),500:50:3500);

subplot(1,2,2);
histogram(E_parameters(:,4),10:10:500);

