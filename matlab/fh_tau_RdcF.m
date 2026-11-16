function [tau, eta_cF] = fh_tau_RdcF(rho_l, fck, dv, fFtud, fFtuk, Dlow, gammaV, fyd)
%FH_TAU_RDCF Annex L (L.19)-(L.20) with EN 1992-1-1:2023 (8.11), (8.16).
if nargin < 6, Dlow = 16; end
if nargin < 7, gammaV = 1.4; end
if nargin < 8, fyd = 435; end
ddg = min(16 + Dlow, 40);
eta_cF = min(max(1.2 - 0.5*fFtuk, 0.4), 1.0);
t = 0.66/gammaV*(100*rho_l*fck*ddg/dv)^(1/3);
tmin = 11/gammaV*sqrt(fck/fyd*ddg/dv);
tau = max(eta_cF*t + fFtud, eta_cF*tmin + fFtud);
end
