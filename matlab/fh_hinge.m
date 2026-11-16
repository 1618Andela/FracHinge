function H = fh_hinge(b, h, law, conc, E, s, rebars, nlay)
%FH_HINGE Builds a cohesive cracked hinge (Olesen 2001) for a rectangular section b x h.
%   s: hinge band width (default h/2); rebars: struct array from fh_rebar (or []).
if nargin < 6 || isempty(s), s = h/2; end
if nargin < 7, rebars = []; end
if nargin < 8, nlay = 400; end
dy = h/nlay;
H.y = ((1:nlay) - 0.5)*dy; H.dA = b*dy*ones(1, nlay);
H.b = b; H.h = h; H.law = law; H.conc = conc; H.E = E; H.s = s; H.rebars = rebars;
H.I = b*h^3/12;
end
