function iSSFP = imag_ssfp_signal(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R,a,b_0,phi_incr)
%IMAG_SSFP_SIGNAL
%    Compute the imaginary component of the bSSFP signal - bSSFP = (M_0r + iM_0i) x (M_x + iM_y)
%    USAGE:
%        ISSFP = IMAG_SSFP_SIGNAL(M_0I,M_0R,RF_DUR,T_1,T_2,T_E,T_R,A,B_0,PHI_INCR)
%
%    M_0I     -> Imaginary component of bSSFP signal (a.u.)
%    M_0I     -> Real component of bSSFP signal (a.u.)
%    RF_DUR   -> Effective duration of applied RF pulse (see SSFP finite RF duration paper)
%    T_1      -> T_1 relaxation time (ms)
%    T_2      -> T_2 relaxation time (ms)
%    T_E      -> Applied Echo time (ms) - typically 0.5*T_R
%    T_R      -> Applied Repetition time (ms)
%    A        -> Applied Flip angle (rad)
%    B_0      -> Induced background field dephasing (rad.T_R) - B_0 = 2*pi*T_R*FieldMap(Hz)
%    PHI_INCR -> Applied rf phase increment at each excitation in order to "move" bSSFP bands %    M_0I     -> Imaginary component of bSSFP signal (a.u.)
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk
%    12-Feb-2017 00:03:27

t2 = 1.0./T_1;
t14 = T_R.*t2;
t3 = exp(-t14);
t4 = 1.0./T_2;
t5 = 1.0./T_R;
t6 = RF_dur.*t5.*(1.0./8.0);
t7 = t6+1.0./8.0;
t8 = T_2.*t2.*t7;
t9 = t8-1.7e1./2.5e1;
t10 = RF_dur.*t9;
t11 = T_R+t10;
t15 = t4.*t11;
t12 = exp(-t15);
t13 = cos(a);
t16 = b_0+phi_incr;
t17 = cos(t16);
t18 = t12.*t17;
t19 = t18-1.0;
t20 = exp(-T_E.*t4);
t21 = sin(a);
t22 = t3-1.0;
t23 = t3.*t13;
t24 = t23-1.0;
t25 = t19.*t24;
t26 = t12-t17;
t27 = t3-t13;
t28 = t25-t12.*t26.*t27;
t29 = 1.0./t28;
iSSFP = M_0r.*t19.*t20.*t21.*t22.*t29-M_0i.*t12.*t20.*t21.*t22.*t29.*sin(t16);
