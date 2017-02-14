function [ costFunction, costFunctionJacobian ] = jsr_cost_function( noisySignal, E, RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr )
%JSR_COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

[ JSRsignal, costFunctionJacobian ] = simulate_JSR_signal( E(2),E(1),RF_dur,E(3),E(4),T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,E(5),phi_incr );

costFunction = JSRsignal(:) - noisySignal(:);


end