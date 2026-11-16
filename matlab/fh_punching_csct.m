function [VR, psi, VRc, VRf] = fh_punching_csct(rc, d, rs, fc, fy, mR, dg, fR1, fR3, Es)
%FH_PUNCHING_CSCT Punching of an interior circular column, Critical Shear Crack Theory
%   (Muttoni 2008) + fibre bridging V_Rf = A_p sigma(psi d / 6) (Maya Duque et al. 2012) with the
%   MC2010 linear cohesive law (fR1, fR3; pass 0, 0 for plain RC). Mean values, N mm MPa.
if nargin < 10, Es = 200000; end
b0 = 2 * pi * (rc + d / 2);
Ap = pi * ((rc + d)^2 - rc^2);
Vflex = 2 * pi * mR * rs / (rs - rc);
sig = @(w) 0;
if fR1 > 0, sig = @(w) fh_sigma_frc(w, fR1, fR3); end
VRc_f = @(p) 0.75 * b0 * d * sqrt(fc) / (1 + 15 * p * d / (16 + dg));
g = @(p) VRc_f(p) + Ap * sig(p * d / 6) - Vflex * (p * d * Es / (1.5 * rs * fy))^(2/3);
phi = 1.5 * rs / d * fy / Es;
if g(phi) > 0
  psi = phi; VR = Vflex;
else
  psi = fzero(g, [1e-9, phi]);
  VR = VRc_f(psi) + Ap * sig(psi * d / 6);
end
VRc = VRc_f(psi); VRf = Ap * sig(psi * d / 6);
end

function s = fh_sigma_frc(w, fR1, fR3)
% MC2010 linear post-cracking law, wu = 2.5 mm, extended linearly to zero
fFts = 0.45 * fR1; fFtu = max(fFts - (fFts - 0.5 * fR3 + 0.2 * fR1), 0);
s = max(fFts - (fFts - fFtu) * w / 2.5, 0);
end
