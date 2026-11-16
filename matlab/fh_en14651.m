function T = fh_en14651(law, conc, E, cmod_max, n)
%FH_EN14651 Simulation of the EN 14651 notched prism (150x150x550, a0 = 25, span 500 mm).
if nargin < 4, cmod_max = 4; end
if nargin < 5, n = 160; end
b = 150; h = 150; a0 = 25; L = 500; hsp = h - a0;
H = fh_hinge(b, hsp, law, conc, E, hsp/2, [], 300);
R = fh_response(H, 2*cmod_max/hsp, n);
I = b*h^3/12;
T.F = 4*R.M/L;
T.CMOD = max(R.theta.*(R.y0 - H.y(1)), 0) + R.theta*a0;
th_el = R.M*H.s/(E*(b*hsp^3/12));
T.delta = T.F*L^3/(48*E*I) + max(R.theta - th_el, 0)*L/4;
cm = [0.5, 1.5, 2.5, 3.5];
for j = 1:4
  Fj = interp1(T.CMOD, T.F, cm(j));
  T.fR(j) = 3*Fj*L/(2*b*hsp^2);
end
mask = T.CMOD <= 0.05; T.FL = max(T.F(mask)); T.fL = 3*T.FL*L/(2*b*hsp^2);
end
