function [N, M, out] = fh_resultants(H, theta, y0)
%FH_RESULTANTS Axial force and bending moment (about mid-depth) for rotation theta and
%   neutral-axis height y0 measured from the tension face.
u = theta*(y0 - H.y);
[sig, w] = fh_layers(H, u);
F = sig.*H.dA; N = sum(F); M = sum(F.*(H.h/2 - H.y));
ss = []; ws = [];
for i = 1:numel(H.rebars)
  rb = H.rebars(i); us = theta*(y0 - rb.y);
  if us <= 0
    sg = max(min(rb.Es*us/H.s, rb.fy), -rb.fy); wsi = 0;
  else
    wsi = interp1(H.y, w, rb.y, 'linear', 'extrap');
    if wsi <= 0
      sg = min(rb.Es*us/H.s, rb.fy);
    else
      sg = fh_rebar_sigma(rb, wsi);
    end
  end
  N = N + rb.As*sg; M = M + rb.As*sg*(H.h/2 - rb.y); ss(end+1) = sg; ws(end+1) = wsi; %#ok<AGROW>
end
out.sig = sig; out.w = w; out.u = u; out.sig_s = ss; out.w_s = ws;
end
