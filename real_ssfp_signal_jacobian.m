function jacobian_rSSFP = real_ssfp_signal_jacobian(M_0i,M_0r,RF_dur,T_1,T_2,T_E,T_R,a,b_0,phi_incr)
%REAL_SSFP_SIGNAL_JACOBIAN
%    Compute the jacobian of the real component of the bSSFP signal - bSSFP = (M_0r + iM_0i) x (M_x + iM_y)
%    USAGE:
%        JACOBIAN_RSSFP = REAL_SSFP_SIGNAL_JACOBIAN(M_0I,M_0R,RF_DUR,T_1,T_2,T_E,T_R,A,B_0,PHI_INCR)
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
%    PHI_INCR -> Applied rf phase increment at each excitation in order to "move" bSSFP bands
%
%    AUTHOR: Rui Pedro A. G. Teixeira - rui.teixeira@kcl.ac.uk
%    12-Feb-2017 00:03:26

t2 = 1.0./T_2;
t3 = 1.0./T_1;
t16 = T_R.*t3;
t4 = exp(-t16);
t5 = 1.0./T_R;
t6 = RF_dur.*t5.*(1.0./8.0);
t7 = t6+1.0./8.0;
t8 = T_2.*t3.*t7;
t9 = t8-1.7e1./2.5e1;
t10 = RF_dur.*t9;
t11 = T_R+t10;
t14 = t2.*t11;
t12 = exp(-t14);
t13 = b_0+phi_incr;
t15 = cos(a);
t17 = cos(t13);
t30 = T_E.*t2;
t18 = exp(-t30);
t19 = sin(a);
t20 = t4-1.0;
t21 = t4.*t15;
t22 = t21-1.0;
t23 = t12.*t17;
t24 = t23-1.0;
t25 = t22.*t24;
t26 = t12-t17;
t27 = t4-t15;
t31 = t12.*t26.*t27;
t28 = t25-t31;
t29 = 1.0./t28;
t32 = 1.0./T_1.^2;
t33 = sin(t13);
t34 = 1.0./t28.^2;
t35 = T_R.*t4.*t12.*t26.*t32;
t40 = t2.*t11.*2.0;
t36 = exp(-t40);
t37 = RF_dur.*t7.*t27.*t32.*t36;
t38 = RF_dur.*t7.*t12.*t26.*t27.*t32;
t39 = t35+t37+t38-T_R.*t4.*t15.*t24.*t32-RF_dur.*t7.*t12.*t17.*t22.*t32;
t41 = 1.0./T_2.^2;
t42 = t11.*t41;
t44 = RF_dur.*t2.*t3.*t7;
t43 = t42-t44;
t45 = t27.*t36.*t43;
t46 = t12.*t26.*t27.*t43;
t47 = t45+t46-t12.*t17.*t22.*t43;
t48 = t12.*t27.*t33;
t49 = t12.*t22.*t33;
t50 = t48+t49;
jacobian_rSSFP = [-t12.*t18.*t19.*t20.*t29.*t33,-t18.*t19.*t20.*t24.*t29,-M_0i.*t18.*t19.*t20.*t24.*t34.*t39-M_0i.*T_R.*t4.*t18.*t19.*t24.*t29.*t32-M_0r.*t12.*t18.*t19.*t20.*t33.*t34.*t39-M_0r.*T_R.*t4.*t12.*t18.*t19.*t29.*t32.*t33-M_0i.*RF_dur.*t7.*t12.*t17.*t18.*t19.*t20.*t29.*t32-M_0r.*RF_dur.*t7.*t12.*t18.*t19.*t20.*t29.*t32.*t33,-M_0i.*t18.*t19.*t20.*t24.*t34.*t47-M_0i.*T_E.*t18.*t19.*t20.*t24.*t29.*t41-M_0i.*t12.*t17.*t18.*t19.*t20.*t29.*t43-M_0r.*t12.*t18.*t19.*t20.*t29.*t33.*t43-M_0r.*t12.*t18.*t19.*t20.*t33.*t34.*t47-M_0r.*T_E.*t12.*t18.*t19.*t20.*t29.*t33.*t41,-M_0r.*t12.*t17.*t18.*t19.*t20.*t29+M_0i.*t12.*t18.*t19.*t20.*t29.*t33-M_0i.*t18.*t19.*t20.*t24.*t34.*t50-M_0r.*t12.*t18.*t19.*t20.*t33.*t34.*t50];
