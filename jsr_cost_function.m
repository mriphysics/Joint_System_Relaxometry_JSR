function [ costFunction, costFunctionJacobian ] = jsr_cost_function( noisySignal, E, RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr )
%USAGE:
% [ costFunction, costFunctionJacobian ] = jsr_cost_function( noisySignal, E, RF_dur,T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,phi_incr );
%
%    JSR_COST_FUNCTION - this function computes the difference between jsr signal model and measured data
%                        it assumes that data is sorted in the correct way [SPGR1;...;SPGRn;rSSFP1;...;rSSFPn;iSSFP1;...;iSSFPn]
%                        The jacobian is also returned in order in case the solver algorithm requires it for performance improvement.
%                        
%    noisySignal -> expected as the concatenation of SPGR and rSSFP and iSSFP signals -  [SPGR1;...;SPGRn;rSSFP1;...;rSSFPn;iSSFP1;...;iSSFPn]
%    E           -> Vector of parameters to be estimated, it assumes:
%                    - E(1) -- rM0, real component of proton density
%                    - E(2) -- iM0, imaginary component of proton density
%                    - E(3) -- T_1 
%                    - E(4) -- T_2
%                    - E(5) -- B0
%   RF_dur       -> Duration of the rf pulse used to produce the excitation (ms)
%   T_E          -> Echo time applied, typically 0.5 T_R_SSFP (ms)
%   T_R_spgr     -> Repetition time of SPGR signal (ms)
%   T_R_ssfp     -> Repetition time of ssfp signal (ms)
%   a_spgr       -> Flip angles used in SPGR acquisition (rad)
%   a_ssfp       -> Flip angles used in SSFP acquisition (rad)
%   phi_incr     -> SSFP rf phase increment that allows to "move" characteristic bands
%
%   Author: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk

[ JSRsignal, costFunctionJacobian ] = simulate_JSR_signal( E(2),E(1),RF_dur,E(3),E(4),T_E,T_R_spgr,T_R_ssfp,a_spgr,a_ssfp,E(5),phi_incr );

costFunction = JSRsignal(:) - noisySignal(:);


end
