% EXAMPLE_HYBRID_BEAM  300x500 beam, 3 phi16 + SFRC class 4c (mean values): M-theta response
fck = 30; fcm = 38; E = 22000*(fcm/10)^0.3; fctm = 0.3*fck^(2/3); GF = 73*fcm^0.18/1000;
fib = fh_law_annexL(4.0, 3.6, fctm, E, 125, 1, 1, 2.5, false);
law = fh_law_envelope(fh_law_mc2010_plain(fctm, GF), fib, 2.5);
conc = fh_concrete('annexL', fcm, E, 4.0);
As = 3*pi*16^2/4;
rb = fh_rebar(40, As, 16, 550, 594, 200000, 0.075, 2*fctm, fctm);
H = fh_hinge(300, 500, law, conc, E, 250, rb, 400);
R = fh_response(H, 0.012, 160);
fprintf('Peak moment (mean laws): %.1f kNm\n', max(R.M)/1e6);
try
  plot(R.theta*1000, R.M/1e6); xlabel('\theta [mrad]'); ylabel('M [kNm]'); grid on;
catch
end
