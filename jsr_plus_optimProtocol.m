%% Obtain Optimized Protocol
clear
Ntrials = 10;
%% Define acquisition limits
TRssfp          = 7;  % as Signal is fairly independent of TR optimization always seeks minimize it.
min_MPRAGE_TR   = 12; %reference FA is 180 sinc-gauss pulse min allows TR is 12ms at 3T for head SAR limits
minTR_spgr      = 7;  %as with ssfp doesn't require 180 reference. min set because of 3.2ms pulse duration


M0 = complex(10,0).*exp(-1i.*(0:0.2:2*pi));

[rM0,iM0,R1,R2,b0,b1] = ndgrid(10,0,[1200],[30],-100:10:100,[0.9:0.1:1.1]);
rf_trms = (1.2*3.2/2.1362);
% idx2delete = find((T1 == 1200 & T2 == 50) |  (T1==800 & T2 == 80));
optimizedFAset  = cell(10, 10, 10, Ntrials);
fval            = cell(10, 10, 10, Ntrials);
pauseFlag = false;

% warning('off','all')
% warning
%%
for Nmprage = [1]
    for Nspgr = [1]
        for Nssfp = [4]
            ParameterSet =[];
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
                if ((Nmprage + Nspgr + Nssfp < 5)||(Nmprage*14 + Nspgr*14 + Nssfp*7 > 57))
                    disp('Less then 6 measurements, ignore iteration');
                    Nmprage_tmp = Nmprage;
                    Nspgr_tmp = Nspgr;
                    Nssfp_tmp = Nssfp;
                    trialN_tmp = trialN;
                    pauseFlag = false;
                    continue;
                end
                
                
                %% For each trial Generate random starting points from Uniform distribution
                TRmprage        = 14;%12 + (15-12).*rand(1,Nmprage);
                alpha_mprage    = deg2rad(6 + (90-2).*rand(1,Nmprage));
                TI              = TRmprage;
                TD              = TRmprage;
                TFEfactor       = 50 + (2000-50).*rand(1,Nmprage);
                t_rage          = TFEfactor.*TRmprage;
                TEmprage        = TRssfp/2;
                
                alpha_ssfp  = deg2rad(6 + (90-6).*rand(1,Nssfp));
                phi_ssfp    = deg2rad(0 + (360-0).*rand(1,Nssfp));
                TEssfp      = TRssfp./2;
                
                alpha_spgr  = deg2rad(6 + (90-2).*rand(1,Nspgr));
                TRspgr      = 7 + (15-7).*rand(1,Nspgr);
                TEspgr      = TRssfp./2;
                
                [optimizedFAset{Nmprage+1,Nspgr+1,Nssfp+1,trialN},fval{Nmprage+1,Nspgr+1,Nssfp+1,trialN}] = JSRplus_FindOptimSet(rM0, iM0, R1, R2, b0, b1,...
                    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
                    TRspgr, TEspgr, alpha_spgr,...
                    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);%M0,T1,T2,b0,rf_phase_incr,TRspgr,alpha_spgr,alpha_ssfp,rf_trms,b1);
                Nmprage_tmp = Nmprage;
                Nspgr_tmp = Nspgr;
                Nssfp_tmp = Nssfp;
                trialN_tmp = trialN;
                pauseFlag = false;
                save('Optimization_BSA_Max_unknownB1')
                ParameterSet = [ParameterSet; optimizedFAset{Nmprage+1,Nspgr+1,Nssfp+1,trialN}'];
                T2 = table(ParameterSet);
                disp(T2);
                disp(ParameterSet)
            end
        end
    end
end
% 
% warning('on','all')
% warning('query','all')
%%
load('Optimization_BSA_Max_unknownB1.mat')
%%
for Nmprage = [0:4]
    for Nspgr = [0:4]
        for Nssfp = [4:9]
            for trialN = 1:Ntrials
                if isempty(fval{Nmprage+1,Nspgr+1,Nssfp+1,trialN})
                    fvalArray(Nmprage+1,Nspgr+1,Nssfp+1,trialN) = NaN;
                else
                    fvalArray(Nmprage+1,Nspgr+1,Nssfp+1,trialN) = fval{Nmprage+1,Nspgr+1,Nssfp+1,trialN};
                end
                
            end
        end
    end
end


for Nmprage = [0:4]
    for Nspgr = [0:4]
        for Nssfp = [4:9]
            [fvalsorted{Nmprage+1,Nspgr+1,Nssfp+1}, sortingTrialIdx{Nmprage+1,Nspgr+1,Nssfp+1}] = sort(squeeze(fvalArray(Nmprage+1,Nspgr+1,Nssfp+1,:)),'ascend');
            for trialN = 1:Ntrials
                rankedProtocol{Nmprage+1,Nspgr+1,Nssfp+1,trialN} = optimizedFAset{Nmprage+1,Nspgr+1,Nssfp+1,sortingTrialIdx{Nmprage+1,Nspgr+1,Nssfp+1}(trialN)};             
            end
        end
    end
end
   
for Nmprage = [0:4]
    for Nspgr = [0:4]
        for Nssfp = [4:9]
            PerformanceMatrix(Nmprage+1,Nspgr+1,Nssfp+1) = fvalsorted{Nmprage+1,Nspgr+1,Nssfp+1}(1);
        end
    end
end
%%
fvalArray = [];
for kk=1:10
    if isempty(fval{2,1,6,kk})
        fvalArray(kk,1) = nan;
    else
        fvalArray(kk,1) = fval{2,1,6,kk};
    end
end
[minTrial,minTrialIdx] = sort(fvalArray,1,'ascend');
optimProtocol = optimizedFAset{2,1,6,minTrialIdx(2)}'