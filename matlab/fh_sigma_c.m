function s = fh_sigma_c(c, eps)
%FH_SIGMA_C Compressive stress (negative) for strain eps (<= 0).
e = -min(eps, 0);
switch c.kind
  case 'elastic-plastic'
    s = -min(c.Ec*e, c.fc); return
  otherwise
    eta = e/c.eps_c1;
    k = 1.05*c.Ec*c.eps_c1/c.fc;
    sig = c.fc*(k*eta - eta.^2)./(1 + (k - 2)*eta);
    sig(eta >= k) = 0;        % branch ends at eta = k (also avoids the pole when k < 2)
    if ~isnan(c.kpost)
      kp = c.kpost;
      sp = c.fc*(kp*eta - eta.^2)./(1 + (kp - 2)*eta);
      sp(eta >= kp) = 0;
      sig(eta > 1) = sp(eta > 1);
    end
    sig(~isfinite(sig)) = 0; sig = min(max(sig, 0), c.fc); sig(e > c.eps_cu) = 0;
    s = -sig;
end
end
