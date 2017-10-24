function [S,J] = simulateJSRplus(M0, R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,varargin)

%% start by generating the MPRAGE signal
% From:
%   Optimization of 3-D MP-RAGE sequences for structural brain imaging. NeuroImage 2000

t0 = abs(M0);
t1 = TRmprage.^-1;
t2 = b1.*alpha_mprage;
t3 = cos(t2);
t4 = log(t3);

R1star  = R1 - t1.*t4;
M0star  = t0.*( (1-exp(-TRmprage.*R1))./(1 - exp(-TRmprage.*R1star)) );

t5 = -t_rage.*R1star;
t6 = exp(t5);

A1      = M0star.*(1 - t6);
B1      = t6;

A2      = t0.*(1 - exp(-TD.*R1));
B2      = exp(-TD.*R1);

t7 = exp(-TI.*R1);
A3      = t0.*(1 - t7);
B3      = cos(pi.*b1).*t7; % Consideres imperfect inversion

A       = A3 + A2.*B3 + A1.*B2.*B3;
B       = B1.*B2.*B3;

M1      = A./(1 - B);

Smprage = (M0star + (M1 - M0star).*exp(t5./2)).*sin(t2).*exp(-TEmprage.*R2);

%% Proceed with SPGR
t8 = exp(-TRspgr.*R1);
t9 = b1.*alpha_spgr;
Sspgr = t0 .* exp(-TEspgr.*R2) .* sin(t9).*((1-t8)./(1 - t8.*cos(t9)));

%% Finalise with bSSFP

t10 = b0.*2.*pi.*10^-3;
t11 = b1.*alpha_ssfp;

b0 = t10.*TRssfp;
b02 = t10.*TEssfp + phi_ssfp;
beta = b0 + phi_ssfp;
cbeta = cos(beta);
sbeta = sin(beta);

calpha = cos(t11);
salpha = sin(t11);

rf_trmsCorrection = 0.68 - 0.125.*(1 + rf_trms./TRssfp).*(R1/R2);

% Obtain acqual steady state solution
E2 = exp(-(TRssfp-rf_trmsCorrection.*rf_trms).*R2);
E1  = exp(-TRssfp.*R1);

d   = (1-E1.*calpha).*(1-E2.*cbeta)-E2.*(E1-calpha).*(E2-cbeta);

Mex = (1-E1).*E2.*salpha.*sbeta./d;
Mey = (1-E1).*((1-E2.*cbeta).*salpha)./d;

% measured SSFP is the transverse magnetization
Sssfp = M0.*exp(1i*b02).*exp(-TEssfp.*R2).*(Mex+1i.*Mey);

%% Now put everything in array to be exported
% S = [real(Smprage);imag(Smprage);real(Sspgr);imag(Sspgr);real(Sssfp);imag(Sssfp)];
S = [(Smprage(:));(Sspgr(:));real(Sssfp(:));imag(Sssfp(:))];

if nargout > 1
    %     J2 = [Jmprage;Jspgr;Jssfp];
    J = simulateJSRplus_Jacobian(S,M0, R1, R2, b0, b1,...
        TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
        TRspgr, TEspgr, alpha_spgr,...
        TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
end


function [J] = simulateJSRplus_Jacobian(S,M0, R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms)

J = zeros(length(S),6);
stepSize    = 1e-3;%sqrt(eps); % differentiation step size
% actual step for each considered parameter
rM0 = real(M0);
iM0 = imag(M0);
if real(M0) == 0
    rM0step = stepSize;
else
    rM0step = rM0.*stepSize;
end

if imag(M0) == 0
    iM0step = stepSize;
else
    iM0step = iM0.*stepSize;
end

R1step      = stepSize;
R2step      = stepSize;

if b0 == 0
    b0step = stepSize;
else
    b0step = b0.*stepSize;
end
b1step       = b1.*stepSize;


J(:,1)  = (simulateJSRplus(rM0 + rM0step + 1i.*iM0, R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./rM0step;

J(:,2)  = (simulateJSRplus(rM0 + 1i.*(iM0+iM0step), R1, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./iM0step;

J(:,3)  = (simulateJSRplus(rM0 + 1i.*iM0, R1+R1step, R2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./R1step;

J(:,4)  = (simulateJSRplus(rM0 + 1i.*iM0, R1, R2+R2step, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./R2step;

J(:,5)  = (simulateJSRplus(rM0 + 1i.*iM0, R1, R2, b0+b0step, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./b0step;

J(:,6)  = (simulateJSRplus(rM0 + 1i.*iM0, R1, R2, b0, b1+b1step,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./b1step;
