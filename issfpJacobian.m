function out1 = issfpJacobian(R1,R2,TEssfp,TRssfp,alpha_ssfp,b0,b1,iM0,phi_ssfp,rM0,rf_trms)
%ISSFPJACOBIAN
%    OUT1 = ISSFPJACOBIAN(R1,R2,TESSFP,TRSSFP,ALPHA_SSFP,B0,B1,IM0,PHI_SSFP,RM0,RF_TRMS)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    20-Sep-2017 12:26:12

t2 = 1.0./R2;
t3 = 1.0./TRssfp;
t4 = rf_trms.*t3.*(1.0./8.0);
t5 = t4+1.0./8.0;
t6 = R1.*t2.*t5;
t7 = t6-1.7e1./2.5e1;
t8 = rf_trms.*t7;
t9 = TRssfp+t8;
t17 = R2.*t9;
t10 = exp(-t17);
t11 = b0+phi_ssfp;
t12 = cos(t11);
t13 = t10.*t12;
t14 = t13-1.0;
t18 = R1.*TRssfp;
t15 = exp(-t18);
t16 = alpha_ssfp.*b1;
t19 = cos(t16);
t32 = R2.*TEssfp;
t20 = exp(-t32);
t21 = sin(t16);
t22 = t15-1.0;
t23 = t15.*t19;
t24 = t23-1.0;
t25 = t14.*t24;
t26 = t10-t12;
t27 = t15-t19;
t33 = t10.*t26.*t27;
t28 = t25-t33;
t29 = 1.0./t28;
t30 = TEssfp.*b0.*pi.*2.0e-3i;
t31 = exp(t30);
t34 = t14.*t20.*t21.*t22.*t29.*1i;
t35 = sin(t11);
t83 = t10.*t20.*t21.*t22.*t29.*t35;
t36 = t34-t83;
t37 = t31.*t36;
t38 = 1.0./t28.^2;
t46 = R2.*t9.*2.0;
t39 = exp(-t46);
t40 = rf_trms.*t5.*t27.*t39;
t41 = TRssfp.*t10.*t15.*t26;
t42 = rf_trms.*t5.*t10.*t26.*t27;
t43 = t40+t41+t42-TRssfp.*t14.*t15.*t19-rf_trms.*t5.*t10.*t12.*t24;
t44 = iM0.*1i;
t45 = rM0+t44;
t48 = R1.*rf_trms.*t2.*t5;
t47 = TRssfp+t8-t48;
t49 = t27.*t39.*t47;
t50 = t10.*t26.*t27.*t47;
t51 = t49+t50-t10.*t12.*t24.*t47;
t52 = TRssfp.^2;
t53 = R2.*t52.*2.0e2;
t54 = rf_trms.^2;
t55 = R1.*t54.*2.5e1;
t56 = R1.*TRssfp.*rf_trms.*2.5e1;
t57 = R2.*TEssfp.*TRssfp.*2.0e2;
t62 = R2.*TRssfp.*rf_trms.*1.36e2;
t58 = t53+t55+t56+t57-t62;
t59 = t3.*t58.*(1.0./2.0e2);
t60 = exp(t59);
t61 = exp(t18);
t63 = t53+t55+t56-t62;
t64 = t3.*t63.*(1.0./2.0e2);
t65 = exp(t64);
t66 = R1.*t52.*2.0e2;
t67 = t53+t55+t56-t62+t66;
t68 = t3.*t67.*(1.0./2.0e2);
t69 = exp(t68);
t70 = R1.*t52.*1.0e2;
t71 = t53+t55+t56-t62+t70;
t72 = t3.*t71.*(1.0./1.0e2);
t73 = exp(t72);
t74 = t3.*t63.*(1.0./1.0e2);
t75 = exp(t74);
t76 = t12.*t65;
t77 = t12.*t19.*t65;
t78 = t12.*t69;
t79 = t12.*t19.*t69;
t80 = t10.*t27.*t35;
t81 = t10.*t24.*t35;
t82 = t80+t81;
out1 = [imag(t37),real(t37),-imag(t31.*t45.*(t14.*t20.*t21.*t22.*t38.*t43.*1i+TRssfp.*t14.*t15.*t20.*t21.*t29.*1i-TRssfp.*t10.*t15.*t20.*t21.*t29.*t35-t10.*t20.*t21.*t22.*t35.*t38.*t43+rf_trms.*t5.*t10.*t12.*t20.*t21.*t22.*t29.*1i-rf_trms.*t5.*t10.*t20.*t21.*t22.*t29.*t35)),-imag(t31.*t45.*(t14.*t20.*t21.*t22.*t38.*t51.*1i+TEssfp.*t14.*t20.*t21.*t22.*t29.*1i-TEssfp.*t10.*t20.*t21.*t22.*t29.*t35+t10.*t12.*t20.*t21.*t22.*t29.*t47.*1i-t10.*t20.*t21.*t22.*t29.*t35.*t47-t10.*t20.*t21.*t22.*t35.*t38.*t51)),imag(-alpha_ssfp.*t45.*exp(TEssfp.*(R2.*1.0e3-b0.*pi.*1i).*(-1.0./5.0e2)).*(t61-1.0).*(exp(t3.*(t53+t55+t56-t62+R2.*TEssfp.*TRssfp.*1.0e2).*(1.0./1.0e2)).*1i-t12.*t60.*1i+t35.*t60).*(t19-t61+t75-t76-t77+t78+t79-t19.*t73).*1.0./(t73+t76+t77-t78-t79+t19.*t61-t19.*t75-1.0).^2),imag(TEssfp.*t31.*t36.*t45.*pi.*2.0e-3i)-imag(t31.*t45.*(t10.*t12.*t20.*t21.*t22.*t29+t10.*t20.*t21.*t22.*t29.*t35.*1i-t14.*t20.*t21.*t22.*t38.*t82.*1i+t10.*t20.*t21.*t22.*t35.*t38.*t82))];
