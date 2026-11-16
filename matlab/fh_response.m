function R = fh_response(H, theta_max, n, Next)
%FH_RESPONSE Moment-rotation response under increasing rotation (monotonic).
if nargin < 3, n = 120; end
if nargin < 4, Next = 0; end
th0 = max(1e-3*H.law.s(1)/H.E*H.s/H.h, 1e-9);
thc = 2*H.law.s(1)*H.s/(H.E*H.h);
th = unique([logspace(log10(th0), log10(theta_max), n), linspace(0.6*thc, 4*thc, max(floor(n/2), 40))]);
th = th(th <= theta_max);
R.theta = [0, th]; R.M = zeros(size(R.theta)); R.y0 = R.M; R.wb = R.M; R.eps_top = R.M;
R.y0(1) = H.h/2; y0 = H.h/2;
for i = 2:numel(R.theta)
  [y0, M, out] = fh_solve(H, R.theta(i), Next, y0);
  R.M(i) = M; R.y0(i) = y0; R.wb(i) = out.w(1); R.eps_top(i) = R.theta(i)*(y0 - H.h)/H.s;
end
end
