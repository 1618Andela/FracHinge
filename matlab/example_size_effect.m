% EXAMPLE_SIZE_EFFECT  Flexural strength of plain concrete vs depth (cracked hinge vs MC2010 5.1-8)
fck = 30; fcm = fck + 8; E = 22000*(fcm/10)^0.3; fctm = 0.3*fck^(2/3); GF = 73*fcm^0.18/1000;
law = fh_law_mc2010_plain(fctm, GF); conc = fh_concrete('ec2', fcm, E);
hs = [100 200 400 800]; ratio = zeros(size(hs)); mc = ratio;
for i = 1:numel(hs)
  H = fh_hinge(100, hs(i), law, conc, E, [], [], 400);
  R = fh_response(H, 40*fctm/E, 120);
  ratio(i) = max(R.M)/(100*hs(i)^2/6)/fctm;
  a = 0.06*hs(i)^0.7/(1 + 0.06*hs(i)^0.7); mc(i) = 1/a;
end
disp([hs; ratio; mc]');
