function law = fh_law_annexL(fR1k, fR3k, fctm, Ecm, lcs, kO, kG, wu, include_matrix)
%FH_LAW_ANNEXL EN 1992-1-1:2023 Annex L law (Figure L.1) written in crack opening w = eps*lcs.
%   Residual branch through (0.5 mm, fFt1,ef) and (2.5 mm, fFt3,ef), (L.7), (L.8).
if nargin < 5 || isempty(lcs), lcs = 125; end
if nargin < 6 || isempty(kO), kO = 1; end
if nargin < 7 || isempty(kG), kG = 1; end
if nargin < 8 || isempty(wu), wu = 2.5; end
if nargin < 9, include_matrix = true; end
f1 = kO*kG*0.37*fR1k; f3 = kO*kG*max(0.57*fR3k - 0.26*fR1k, 0);
line = @(w) max(f1 + (f3 - f1)*(w - 0.5)/2, 0);
w0 = fctm/Ecm*lcs;
if include_matrix
  law.w = [0, w0, wu]; law.s = [fctm, line(w0), line(wu)];
else
  law.w = [0, wu]; law.s = [line(0), line(wu)];
end
law.tail = 'zero'; law.name = 'EC2 Annex L';
end
