function [sig, w] = fh_layers(H, u)
%FH_LAYERS Stresses and crack openings for band elongations u (row vector), vectorised.
%   Cracked layer: u = s*sigma/E + w with sigma = sigma_w(w) (fictitious crack, Olesen 2001).
s = H.s; E = H.E; law = H.law;
sig = zeros(size(u)); w = zeros(size(u));
comp = u <= 0;
sig(comp) = fh_sigma_c(H.conc, u(comp)/s);
ucr = s*law.s(1)/E;
el = (~comp) & (u <= ucr);
sig(el) = E*u(el)/s;
cr = (~comp) & (u > ucr);
if any(cr)
  uc = u(cr)';                                  % column
  w0 = law.w(1:end-1); w1 = law.w(2:end);
  bb = (law.s(2:end) - law.s(1:end-1))./(w1 - w0); aa = law.s(1:end-1) - bb.*w0;
  den = 1 + s*bb/E;                             % row
  WW = (uc - s*aa/E)./den;                      % ncr x nseg (implicit expansion)
  ok = (den > 0) & (WW >= w0 - 1e-12) & (WW <= w1 + 1e-12);
  WW(~ok) = Inf;
  wm = min(WW, [], 2);
  if strcmp(law.tail, 'zero'), wt = uc; else, wt = uc - s*law.s(end)/E; end
  wm(~isfinite(wm)) = wt(~isfinite(wm));
  w(cr) = wm'; sig(cr) = fh_sigma_w(law, wm');
end
end
