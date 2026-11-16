function [Vu, xd, xi] = fh_masm_sfrc(b, d, a, rho, fcm, sigma_tu, Es)
%FH_MASM_SFRC Shear strength of SFRC beams without stirrups, Multi-Action Shear Model
%   (Mari et al. 2020, Eq. 26), mean values. Units N, mm, MPa.
%   sigma_tu: residual tensile stress of the cohesive law (e.g. fFtu).
if nargin < 7, Es = 200000; end
fct = 0.3 * min(fcm - 8, 60)^(2/3);
Ec = 22000 * (fcm / 10)^0.3;
nr = Es / Ec * rho;
xd = nr * (-1 + sqrt(1 + 2 / nr));
xi = max(2 / sqrt(1 + max(d, 100) / 200) * (d / a)^0.2, 0.45);
s = min(sigma_tu / fct, 1);
Vu = (xi * (0.84 - 0.10 * s) * xd + 0.08 + 1.10 * s) * fct * b * d;
end
