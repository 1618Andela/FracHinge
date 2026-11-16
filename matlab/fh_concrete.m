function c = fh_concrete(kind, fcm, Ecm, fR1k)
%FH_CONCRETE Compression laws: 'ec2' (Formulae 5.6-5.10), 'annexL' ((L.12)-(L.14)),
%   'elastic-plastic'. Returns a struct used by fh_sigma_c.
if nargin < 3 || isempty(Ecm), Ecm = 22000*(fcm/10)^0.3; end
c.kind = kind; c.fc = fcm; c.Ec = Ecm; c.kpost = NaN;
switch kind
  case 'ec2'
    c.eps_c1 = min(0.7*fcm^(1/3), 2.8)/1000;
    c.eps_cu = min(2.8 + 14*(1 - fcm/108)^4, 3.5)/1000;
  case 'annexL'
    c.eps_c1 = 0.7*fcm^(1/3)*(1 + 0.03*fR1k)/1000;
    c.kpost = 1 + 20/sqrt(82 - 2.2*fR1k);
    c.eps_cu = c.kpost*c.eps_c1;
  case 'elastic-plastic'
    c.eps_c1 = fcm/Ecm; c.eps_cu = Inf;
end
end
