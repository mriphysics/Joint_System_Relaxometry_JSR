function [ JSRsignal, JSRsignalJacobian ] = simulate_JSR_signal( M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,b_0,phi_incr )
%SIMULATE_JSR_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

JSRsignal = [spgr_signal(M_0i,M_0r,T_1,T_2,T_E,T_R_spgr,a_spgr(:));...
    real_ssfp_signal(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_ssfp,a_ssfp(:),b_0(:),phi_incr(:));...
    imag_ssfp_signal(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_ssfp,a_ssfp(:),b_0(:),phi_incr(:))];

if nargout > 1   % Two output arguments
   JSRsignalJacobian = [spgr_signal_jacobian(M_0i,M_0r,T_1,T_2,T_E,T_R_spgr,a_spgr(:));...
       real_ssfp_signal_jacobian(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_ssfp,a_ssfp(:),b_0(:),phi_incr(:));...
       imag_ssfp_signal_jacobian(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R_ssfp,a_ssfp(:),b_0(:),phi_incr(:))];   % Jacobian of the function evaluated at x
end

end