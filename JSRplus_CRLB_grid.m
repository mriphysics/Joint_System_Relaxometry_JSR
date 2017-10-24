function [CF] = JSRplus_CRLB_grid(rM0grid, iM0grid, R1grid, R2grid, b0grid, b1grid,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms)

noise_var = 0.002;
covM = noise_var*eye(length([alpha_mprage(:);alpha_spgr(:);alpha_ssfp(:);alpha_ssfp(:)]));
gridL = length(rM0grid(:));
CF = zeros(gridL,1);

parfor ii=1:gridL
    
    [~,J] = simulateJSRplus(rM0grid(ii) + 1i.*iM0grid(ii), R1grid(ii), R2grid(ii), b0grid(ii), b1grid(ii),...
        TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
        TRspgr, TEspgr, alpha_spgr,...
        TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
    
    F = J'/covM*J;
    CRLB = (diag(F\eye(size(F))));
    
   CF(ii) = sqrt(sum(CRLB(3:4)./[R1grid(ii).^2;R2grid(ii).^2]));
   
end
CF = max(CF);