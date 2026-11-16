function [Vu, xd, xi, lam, eta] = fh_masm_sfrc(b, d, a, rho, fcm, sigma_tu, Es, size_effect, aligned)
%FH_MASM_SFRC Shear strength of SFRC beams without stirrups, Multi-Action Shear Model
%   (Mari et al. 2020, Eq. 26), mean values. Units N, mm, MPa.
%   sigma_tu: residual tensile stress of the cohesive law (e.g. fFtu).
%   size_effect (default true): fibre term times Bazant/ACI 318-19 lambda_s = sqrt(2/(1+d/254)) <= 1.
if nargin < 7 || isempty(Es), Es = 200000; end
if nargin < 8 || isempty(size_effect), size_effect = true; end
if nargin < 9, aligned = true; end
%   aligned (default true): fibres aligned with the beam axis, bridging stress x sin(theta_cr),
%   cot(theta_cr) = 0.85 d/(d - x) <= 2.5 (FracHinge 1.3); false = isotropic orientation.
fct = 0.3 * min(fcm - 8, 60)^(2/3);
Ec = 22000 * (fcm / 10)^0.3;
nr = Es / Ec * rho;
xd = nr * (-1 + sqrt(1 + 2 / nr));
xi = max(2 / sqrt(1 + max(d, 100) / 200) * (d / a)^0.2, 0.45);
cot_t = min(0.85 / (1 - xd), 2.5);
eta = 1;
if aligned, eta = 1 / sqrt(1 + cot_t^2); end
s = min(eta * sigma_tu / fct, 1);
lam = 1;
if size_effect, lam = min(sqrt(2 / (1 + d / 254)), 1); end
Vu = (xi * (0.84 - 0.10 * s) * xd + 0.08 + 1.10 * s * lam) * fct * b * d;
end
