function [CF] = JSRplus_CRLB_grid(rM0grid, iM0grid, T1grid, T2grid, b0grid, b1grid,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,knownB1)
noise_var = 0.0008;
covM = noise_var*eye(length([alpha_mprage(:);alpha_spgr(:);alpha_ssfp(:);alpha_ssfp(:)]));
gridL = length(rM0grid(:));
CF = zeros(gridL,1);

if knownB1
    parfor ii=1:gridL
        [~,J] = simulateJSRplus(rM0grid(ii) + 1i.*iM0grid(ii), T1grid(ii), T2grid(ii), b0grid(ii), b1grid(ii),...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        
        J(:,6) = [];
        F = J'/covM*J;
        CRLB = (diag(F\eye(size(F))));
        
        CF(ii) = sqrt(sum(CRLB(3:4)./[T1grid(ii).^2;T2grid(ii).^2]));
    end
else
    
    parfor ii=1:gridL
        
        [~,J] = simulateJSRplus(rM0grid(ii) + 1i.*iM0grid(ii), T1grid(ii), T2grid(ii), b0grid(ii), b1grid(ii),...
            TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
            TRspgr, TEspgr, alpha_spgr,...
            TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
        
        F = J'/covM*J;
        CRLB = (diag(F\eye(size(F))));
        
        CF(ii) = sqrt(sum(CRLB(3:4)./[T1grid(ii).^2;T2grid(ii).^2]));
        
    end
end
CF = max(CF);