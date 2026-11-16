function [y0, M, out] = fh_solve(H, theta, Next, y0g)
%FH_SOLVE Neutral axis satisfying N(theta, y0) = Next (continuation from y0g).
if nargin < 3, Next = 0; end
if nargin < 4 || isempty(y0g), y0g = H.h/2; end
f = @(y) fh_resultants(H, theta, y) - Next;
fr = f(y0g); y0 = NaN;
if fr == 0, y0 = y0g; end
step = 0.01*H.h;
while isnan(y0) && step < 2*H.h
  lo = y0g - step; hi = y0g + step;
  if f(lo)*fr < 0, y0 = fzero(f, [lo, y0g]); break; end
  if f(hi)*fr < 0, y0 = fzero(f, [y0g, hi]); break; end
  step = 1.8*step;
end
if isnan(y0)
  g = linspace(-0.5*H.h, 1.5*H.h, 161); v = arrayfun(f, g);
  k = find(v(1:end-1).*v(2:end) < 0, 1);
  if isempty(k), error('fh_solve: no equilibrium at theta = %g', theta); end
  y0 = fzero(f, [g(k), g(k+1)]);
end
[~, M, out] = fh_resultants(H, theta, y0);
end
