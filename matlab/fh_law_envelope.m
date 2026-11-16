function env = fh_law_envelope(lawA, lawB, w_max)
%FH_LAW_ENVELOPE Upper envelope max(sigma_A, sigma_B) (MC2010 5.6.5 rule).
if nargin < 3, w_max = max(lawA.w(end), lawB.w(end)); end
w = unique([lawA.w, lawB.w, logspace(-4, log10(w_max), 400)]);
w = w(w <= w_max); if w(1) > 0, w = [0, w]; end
s = max(fh_sigma_w(lawA, w), fh_sigma_w(lawB, w));
env.w = w; env.s = s; env.tail = 'zero'; env.name = 'envelope';
end
