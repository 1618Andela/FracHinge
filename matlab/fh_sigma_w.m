function s = fh_sigma_w(law, w)
%FH_SIGMA_W Evaluates a piecewise-linear cohesive law sigma(w).
s = interp1(law.w, law.s, min(max(w, 0), law.w(end)), 'linear');
if strcmp(law.tail, 'zero')
  s(w > law.w(end)) = 0;
else
  s(w > law.w(end)) = law.s(end);
end
end
