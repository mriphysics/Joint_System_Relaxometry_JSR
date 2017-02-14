function jacobian = spgr_signal_jacobian(M_0i,M_0r,T_1,T_2,T_E,T_R,a)
%SPGR_SIGNAL_JACOBIAN
%    JACOBIAN = SPGR_SIGNAL_JACOBIAN(M_0I,M_0R,T_1,T_2,T_E,T_R,A)
%
%    M_0i -> Imaginary component of proton density
%    M_0r -> Real component of proton density
%    T_1  -> T_1 relaxation time (ms)
%    T_2  -> T_2 relaxation time (ms)
%    T_E  -> Applied Echo Time (ms)
%    T_R  -> Repetition time between excitations (ms)
%    a    -> Flip angle of excitation (rad)
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk
%    12-Feb-2017 00:04:46

t2 = 1.0./T_1;
t7 = T_R.*t2;
t3 = exp(-t7);
t4 = 1.0./T_2;
t17 = T_E.*t4;
t5 = exp(-t17);
t6 = sin(a);
t8 = t3-1.0;
t9 = cos(a);
t10 = t3.*t9;
t11 = t10-1.0;
t12 = 1.0./t11;
t13 = M_0i.^2;
t14 = M_0r.^2;
t15 = t13+t14;
t16 = 1.0./sqrt(t15);
t18 = 1.0./T_1.^2;
t19 = sqrt(t15);
jacobian = [M_0r.*t5.*t6.*t8.*t12.*t16,M_0i.*t5.*t6.*t8.*t12.*t16,T_R.*t3.*t5.*t6.*t12.*t18.*t19-T_R.*t3.*t5.*t6.*t8.*t9.*1.0./t11.^2.*t18.*t19,1.0./T_2.^2.*T_E.*t5.*t6.*t8.*t12.*t19,zeros(length(a),1)];
