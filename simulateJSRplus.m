function [S,J] = simulateJSRplus(M0, T1, T2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms,varargin)

%% start by generating the MPRAGE signal
% From:
%   Optimization of 3-D MP-RAGE sequences for structural brain imaging. NeuroImage 2000
if ~isempty(alpha_mprage)
    
    t0 = abs(M0);
    t1 = TRmprage.^-1;
    t2 = b1.*alpha_mprage;
    t3 = cos(t2);
    t4 = log(t3);
    
    R1star  = 1./T1 - t1.*t4;
    M0star  = t0.*( (1-exp(-TRmprage./T1))./(1 - exp(-TRmprage.*R1star)) );
    
    t5 = -t_rage.*R1star;
    t6 = exp(t5);
    
    A1      = M0star.*(1 - t6);
    B1      = t6;
    
    A2      = t0.*(1 - exp(-TD./T1));
    B2      = exp(-TD./T1);
    
    t7 = exp(-TI./T1);
    A3      = t0.*(1 - t7);
    B3      = cos(pi.*b1).*t7; % Consideres imperfect inversion
    
    A       = A3 + A2.*B3 + A1.*B2.*B3;
    B       = B1.*B2.*B3;
    
    M1      = A./(1 - B);
    
    Smprage = (M0star + (M1 - M0star).*exp(t5./2)).*sin(t2).*exp(-TEmprage./T2);
else
    Smprage = [];
end
%% Proceed with SPGR
if ~isempty(alpha_spgr)
    E1 = exp(-TRspgr(:)./T1);
    % exp(-TE./T2') is assumed in M0, this requires same echo time to be used
    % between ssfp and spgr acquisition;
    Sspgr = ((abs(M0).*exp(-TEspgr./T2)).*sin(b1.*alpha_spgr(:)).*(1-E1(:)))./(1-E1(:).*cos(b1.*alpha_spgr(:)));
else
    Sspgr = [];
end
%% Finalise with bSSFP
b02 = b0.*2.*pi.*TEssfp.*10^-3 ;
b0 = b0.*2.*pi.*TRssfp.*10^-3; % convert from Hz to radians of precession per TR
% amount of precession due to b0 inhomogeniety and increment of rf applied
% rf phase results in a some of (angles) page 212 Magnetic Reonance
% Imaging - Theory and Practice springer, M.T. Vlaardingerbroek and J.A.
% den Boer;
beta = b0+phi_ssfp(:);
cbeta = cos(beta(:));
sbeta = sin(beta(:));

% prcompute cos(alpha) and sin(alpha) for each one of the measurements
calpha = cos(b1.*alpha_ssfp(:));
salpha = sin(b1.*alpha_ssfp(:));

% SSFP signal with finite RF pulses, O. Biere and K. Scheffler, 2009
% assumes spins assume as less T2 decay then traditionally expected as
% there is no T2 decay when relaxation is "going through Z"
rf_trmsCorrection = 0.68 - 0.125.*(1 + rf_trms(:)./TRssfp(:)).*(T2./T1);


% Obtain acqual steady state solution
% E2  = exp(-TRssfp(:).*eR2(:));
E2 = exp(-(TRssfp-rf_trmsCorrection.*rf_trms(:))./T2);
E1  = exp(-TRssfp(:)./T1);

d   = (1-E1.*calpha).*(1-E2.*cbeta)-E2.*(E1-calpha).*(E2-cbeta);

Me = [(1-E1).*E2.*salpha.*sbeta./d (1-E1).*((1-E2.*cbeta).*salpha)./d].*exp(-TEssfp./T2);

% measured SSFP is the transverse magnetization
Sssfp = M0*exp(-1i.*b02).*[Me(:,1)+1i.*Me(:,2)];

%% Now put everything in array to be exported
% S = [real(Smprage);imag(Smprage);real(Sspgr);imag(Sspgr);real(Sssfp);imag(Sssfp)];
S = [(Smprage(:));(Sspgr(:));real(Sssfp(:));imag(Sssfp(:))];

if nargout > 1
    %     J2 = [Jmprage;Jspgr;Jssfp];
    J = simulateJSRplus_Jacobian(S,M0, T1, T2, b0, b1,...
        TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
        TRspgr, TEspgr, alpha_spgr,...
        TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms);
end


function [J] = simulateJSRplus_Jacobian(S,M0, T1, T2, b0, b1,...
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
    rM0step = abs(stepSize);
end

if imag(M0) == 0
    iM0step = stepSize;
else
    iM0step = abs(stepSize);
end

T1step      = stepSize;
T2step      = stepSize;

if b0 == 0
    b0step = stepSize;
else
    b0step = abs(stepSize);
end
b1step       = stepSize;


J(:,1)  = (simulateJSRplus(rM0 + rM0step + 1i.*iM0, T1, T2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./rM0step;

J(:,2)  = (simulateJSRplus(rM0 + 1i.*(iM0+iM0step), T1, T2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./iM0step;

J(:,3)  = (simulateJSRplus(rM0 + 1i.*iM0, T1+T1step, T2, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./T1step;

J(:,4)  = (simulateJSRplus(rM0 + 1i.*iM0, T1, T2+T2step, b0, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./T2step;

J(:,5)  = (simulateJSRplus(rM0 + 1i.*iM0, T1, T2, b0+b0step, b1,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./b0step;

J(:,6)  = (simulateJSRplus(rM0 + 1i.*iM0, T1, T2, b0, b1+b1step,...
    TRmprage, TEmprage, alpha_mprage, t_rage, TI, TD, ...
    TRspgr, TEspgr, alpha_spgr,...
    TRssfp, TEssfp, alpha_ssfp, phi_ssfp, rf_trms) - S)./b1step;
