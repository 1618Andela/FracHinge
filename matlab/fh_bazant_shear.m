function [Vn, Vc, VF, Vs, c, gd] = fh_bazant_shear(b, d, rho_l, fc, sigma_tu, rho_w, fyw, stirrups_min, aligned, two_way, Es)
%FH_BAZANT_SHEAR Shear (one-way) or punching (two-way) strength, unified approach of Frosch, Yu,
%   Cusatis & Bazant (2017) with the energetic size-effect law (ACI 446), plus cohesive fibre bridging.
%   Units N, mm, MPa. b = web width (one-way) or control perimeter b0 at d/2 (two_way = true).
%   V_c = 5 sqrt(f'c[psi]) b c gamma_d (one-way), 10 sqrt(f'c[psi]) b0 c gamma_d (two-way),
%   c = k d (elastic cracked neutral axis), gamma_d = 1.4/sqrt(1+d/d0) <= 1 (d0 = 254 / 2540 mm),
%   V_F = eta lambda_s sigma_tu b (d - c), lambda_s = sqrt(2/(1+d/254)) <= 1, eta = 1/sqrt(2) (aligned),
%   V_s = rho_w b fyw d, V_n = min(V_c + V_F + V_s, 5 (V_c + V_F)) (one-way).
if nargin < 5 || isempty(sigma_tu), sigma_tu = 0; end
if nargin < 6 || isempty(rho_w), rho_w = 0; end
if nargin < 7 || isempty(fyw), fyw = 500; end
if nargin < 8 || isempty(stirrups_min), stirrups_min = false; end
if nargin < 9 || isempty(aligned), aligned = true; end
if nargin < 10 || isempty(two_way), two_way = false; end
if nargin < 11 || isempty(Es), Es = 200000; end
rn = rho_l * Es / (4700 * sqrt(fc));
c = (sqrt(2*rn + rn^2) - rn) * d;
if stirrups_min, d0 = 2540; else, d0 = 254; end
if d < d0, gd = 1; else, gd = min(1.4 / sqrt(1 + d/d0), 1); end
kv = 5; if two_way, kv = 10; end
Vc = kv / 12.0416 * sqrt(fc) * b * c * gd;
lam = min(sqrt(2 / (1 + d/254)), 1);
if aligned, eta = 1/sqrt(2); else, eta = 1; end
VF = eta * lam * sigma_tu * b * max(d - c, 0);
if two_way
  Vs = 0; Vn = Vc + VF;
else
  Vs = rho_w * b * fyw * d; Vn = min(Vc + VF + Vs, 5 * (Vc + VF));
end
end
