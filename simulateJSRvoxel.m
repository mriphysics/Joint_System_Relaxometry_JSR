%% Simulate the signals necessary for the JSR fit
%

function [signal] = simulateJSRvoxel(M0,T1,T2,b0,rf_phase_incr,TRspgr,TRssfp,TE,alpha_spgr,alpha_ssfp,rf_trms,ssfpOnly)
if T2>T1
    signal = 1e15*ones(length([alpha_spgr(:);alpha_ssfp(:)]),1);
    return;
end
if nargin == 0
    
    fprintf('Test mode \n');
    fprintf('Creating standard variables\n');
    M0 = 1;
    T1 = 850;
    T2 = 40;
    b0 = 0.4;
    rf_phase_incr = zeros(2*90,1);
    rf_phase_incr(1:90) = pi;
    TRspgr = 8.5;
    TRssfp = 4;
    TE = TRssfp/2;
    alpha_spgr = 1:20;
    alpha_ssfp = [1:90 1:90];
    rf_trms = 1;
    
end
if nargin < 12
    ssfpOnly = false;
end

%% Compute SPGR signal equations
% allow to have different TR for different measurements
% If only one TRspgr is used then assume same TRspgr was required for all
% FA measurements;

Nspgr = length(alpha_spgr);
if length(TRspgr) == 1
    TRspgr = TRspgr.*ones(Nspgr,1);
end
E1 = exp(-TRspgr(:)./T1);

% exp(-TE./T2) is assumed in M0, this requires same echo time to be used
% between ssfp and spgr acquisition;
SPGRsignal = ((abs(M0).*exp(-TE./T2)).*sind(alpha_spgr(:)).*(1-E1(:)))./(1-E1(:).*cosd(alpha_spgr(:)));

%% Compute bSSFP signal equations
% allow to have different TR for different measurements
% If only one TRssfp is used then assume same TRssfp was required for all
% FA measurements;

Nssfp = length(alpha_ssfp);
if length(TRssfp) == 1
    TRssfp = TRssfp.*ones(Nssfp,1);
end
if length(rf_phase_incr) == 1
    rf_phase_incr = rf_phase_incr.*ones(Nssfp,1);
end
if length(rf_trms) == 1;
    rf_trms = rf_trms.*ones(Nssfp,1);
end
if length(b0) == 1
    b0 = b0.*ones(Nssfp,1);
end

b0 = b0.*2.*pi.*TRssfp.*10^-3; % convert from Hz to radians of precession per TR
% amount of precession due to b0 inhomogeniety and increment of rf applied
% rf phase results in a some of (angles) page 212 Magnetic Reonance
% Imaging - Theory and Practice springer, M.T. Vlaardingerbroek and J.A.
% den Boer;
beta = b0+rf_phase_incr(:);
cbeta = cos(beta(:));
sbeta = sin(beta(:));

% prcompute cos(alpha) and sin(alpha) for each one of the measurements
calpha = cosd(alpha_ssfp(:));
salpha = sind(alpha_ssfp(:));

% SSFP signal with finite RF pulses, O. Biere and K. Scheffler, 2009
% assumes spins assume as less T2 decay then traditionally expected as
% there is no T2 decay when relaxation is "going through Z"
rf_trmsCorrection = 0.68 - 0.125.*(1 + rf_trms(:)./TRssfp(:)).*(T2./T1);

% Obtain acqual steady state solution
E2 = exp(-(TRssfp-rf_trmsCorrection.*rf_trms(:))./T2);
E1  = exp(-TRssfp(:)./T1);

d   = (1-E1.*calpha).*(1-E2.*cbeta)-E2.*(E1-calpha).*(E2-cbeta);

Me = [(1-E1).*E2.*salpha.*sbeta./d (1-E1).*((1-E2.*cbeta).*salpha)./d].*exp(-TE./T2);

% measured SSFP is the transverse magnetization
SSFPsignal = M0*[Me(:,1)+1i.*Me(:,2)];

signal = [(SPGRsignal(:));SSFPsignal(:)];
if ssfpOnly
    signal = [SSFPsignal(:)];
end

end
