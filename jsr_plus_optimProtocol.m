%% Obtain Optimized Protocol
clear
Ntrials = 1;
%% Define acquisition limits
TRssfp          = 7;  % as Signal is fairly independent of TR optimization always seeks minimize it.
min_MPRAGE_TR   = 12; %reference FA is 180 sinc-gauss pulse min allows TR is 12ms at 3T for head SAR limits
minTR_spgr      = 7;  %as with ssfp doesn't require 180 reference. min set because of 3.2ms pulse duration


M0 = complex(10,0);
[rM0,iM0,R1,R2,b0,b1] = ndgrid(1e1,1e1,[1./2200],[1./180],-150:1:150,[0.8:0.1:1.2]);
rf_trms = (1.2*3.2/2.1362);
% idx2delete = find((T1 == 1200 & T2 == 50) | (T1==800 & T2 == 80));
% optimizedFAset  = cell(10, 10, 10, Ntrials);
% fval            = cell(10, 10, 10, Ntrials);
pauseFlag = false;

% warning('off','all')
% warning
%%
for Nmprage = [0:9]
    for Nspgr = [0:9]
        for Nssfp = [2:9]
            for trialN = 1:Ntrials
                if pauseFlag
                    Nmprage = Nmprage_tmp;
                    Nspgr = Nspgr_tmp;
                    Nssfp = Nssfp_tmp;
                    trialN = trialN_tmp;
                    pauseFlag = false;
                end
                pauseFlag = true; % after testing if stopped assume something went wrong
                
                T = table(Nmprage,Nspgr,Nssfp,trialN);
                disp(T);
                if ~isempty(optimizedFAset{Nspgr+1,Nssfp+1,trialN})
                    continue;
                end
                if Nmprage + Nspgr + Nssfp < 6
                    disp('Less then 6 measurements, ignore iteration');
                    Nmprage_tmp = Nmprage;
                    Nspgr_tmp = Nspgr;
                    Nssfp_tmp = Nssfp;
                    trialN_tmp = trialN;
                    pauseFlag = false;
                    continue;
                end
                
                %% For each trial Generate random starting points from Uniform distribution
                TRmprage        = 12 + (100-12).*rand(1,Nmprage);
                alpha_mprage    = deg2rad(6 + (90-2).*rand(1,Nmprage));
                TI              = TRmprage;
                TD              = TRmprage;
                TFEfactor       = 50 + (1000-50).*rand(1,Nmprage);
                t_rage          = TFEfactor.*TRmprage;
                TEmprage        = TRssfp/2;
                
                alpha_ssfp  = deg2rad(6 + (90-6).*rand(1,Nssfp));
                phi_ssfp    = 0 + (2*pi - 0).*rand(1,Nssfp);
                TEssfp      = TRssfp./2;
                
                alpha_spgr  = deg2rad(6 + (90-2).*rand(1,Nspgr));
                TRspgr      = 12 + (100-12).*rand(1,Nspgr);
                TEspgr      = TRssfp./2;
                
                [optimizedFAset{Nmprage+1,Nspgr+1,Nssfp+1,trialN},fval{Nspgr+1,Nssfp+1,trialN}] = JSRplus_FindOptimSet(rM0, iM0, R1, R2, b0, b1,...
                    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
                    TRspgr, TEspgr, alpha_spgr,...
                    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);%M0,T1,T2,b0,rf_phase_incr,TRspgr,alpha_spgr,alpha_ssfp,rf_trms,b1);
                
                Nmprage_tmp = Nmprage;
                Nspgr_tmp = Nspgr;
                Nssfp_tmp = Nssfp;
                trialN_tmp = trialN;
                pauseFlag = false;
                save('Optimization')
            end
        end
    end
end
% 
% warning('on','all')
% warning('query','all')