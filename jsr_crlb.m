function [crlb] = jsr_crlb(M_0i, M_0r, RF_dur, T_1, T_2, T_E, T_R_spgr, T_R_ssfp, a_spgr, a_ssfp, b_0, phi_incr, noise_std)

% Compute the jacobian
[~,jsr_jacobian] = simulate_JSR_signal(M_0i, M_0r, RF_dur, T_1, T_2, T_E, T_R_spgr, T_R_ssfp, a_spgr, a_ssfp, b_0, phi_incr);

% Compute Fisher information matrix
FisherM = jsr_jacobian.'/diag(noise_std*ones(length([a_spgr(:);a_ssfp(:);a_ssfp]),1),0)*jsr_jacobian; % a_ssfp is repeated because we make use of both real and imaginary component of ssfp signal

% the CRLB is the inverse of the Fisher matrix
crlb = diag(inv(FisherM)); % we are usually only interested in the diagonal elements of crlb
end
