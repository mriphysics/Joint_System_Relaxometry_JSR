function spgr = spgr_signal(M_0i,M_0r,T_1,T_2,T_E,T_R,a)
%SPGR_SIGNAL
%    SPGR = SPGR_SIGNAL(M_0I,M_0R,T_1,T_2,T_E,T_R,A)
%
%    M_0i -> Imaginary component of proton density
%    M_0r -> Real component of proton density
%    T_1  -> T_1 relaxation time (ms)
%    T_2  -> T_2 relaxation time (ms)
%    T_E  -> Applied Echo Time (ms)
%    T_R  -> Repetition time between excitations (ms)
%    a    -> Flip angle of excitation (rad)
%
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk
%    12-Feb-2017 00:04:45

t2 = 1.0./T_1;
t3 = exp(-T_R.*t2);
spgr = (exp(-T_E./T_2).*sin(a).*sqrt(M_0i.^2+M_0r.^2).*(t3-1.0))./(t3.*cos(a)-1.0);
