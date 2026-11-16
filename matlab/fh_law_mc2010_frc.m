function law = fh_law_mc2010_frc(fR1, fR3, wu, w_end)
%FH_LAW_MC2010_FRC MC2010 linear post-cracking law, Eqs. (5.6-5), (5.6-6).
%   fFts = 0.45 fR1; fFtu = fFts - wu/CMOD3 (fFts - 0.5 fR3 + 0.2 fR1) >= 0, CMOD3 = 2.5 mm.
if nargin < 3, wu = 2.5; end
fFts = 0.45*fR1;
fFtu = max(fFts - wu/2.5*(fFts - 0.5*fR3 + 0.2*fR1), 0);
if nargin < 4 || w_end <= wu
  law.w = [0, wu]; law.s = [fFts, fFtu];
else
  slope = (fFtu - fFts)/wu;
  law.w = [0, wu, w_end]; law.s = [fFts, fFtu, max(fFts + slope*w_end, 0)];
end
law.tail = 'zero'; law.name = 'MC2010 FRC linear';
end
