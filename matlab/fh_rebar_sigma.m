function s = fh_rebar_sigma(rb, w)
%FH_REBAR_SIGMA Bar stress at the crack for a crack opening w (pull-out on both faces).
w = max(w, 0);
wy = rb.fy^2*rb.phi/(4*rb.tau0*rb.Es);
if w <= wy
  s = sqrt(4*rb.tau0*rb.Es*w/rb.phi);
else
  c = (w - wy)*4*rb.tau1/(2*rb.phi); A = 1/(2*rb.Esh); B = rb.fy/rb.Es;
  s = rb.fy + (-B + sqrt(B^2 + 4*A*c))/(2*A);
end
s = min(s, rb.fu);
end
