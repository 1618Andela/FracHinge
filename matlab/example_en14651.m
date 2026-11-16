% EXAMPLE_EN14651  Forward simulation of the EN 14651 test with code laws (class 4c)
fck = 30; fcm = 38; E = 22000*(fcm/10)^0.3; fctm = 0.3*fck^(2/3); GF = 73*fcm^0.18/1000;
conc = fh_concrete('elastic-plastic', fcm, E);
lawA = fh_law_envelope(fh_law_mc2010_plain(fctm, GF), fh_law_mc2010_frc(4.0, 3.6, 2.5, 4.5));
lawB = fh_law_annexL(4.0, 3.6, fctm, E, 125, 1, 1, 4.5, true);
TA = fh_en14651(lawA, conc, E, 4, 140); TB = fh_en14651(lawB, conc, E, 4, 140);
fprintf('MC2010 law : fR = %s  fL = %.2f\n', mat2str(TA.fR, 3), TA.fL);
fprintf('Annex L law: fR = %s  fL = %.2f\n', mat2str(TB.fR, 3), TB.fL);
try
  plot(TA.CMOD, TA.F/1000, TB.CMOD, TB.F/1000); xlabel('CMOD [mm]'); ylabel('F [kN]');
  legend('MC2010 linear + plain', 'EC2 Annex L'); grid on;
catch
end
