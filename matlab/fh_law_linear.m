function law = fh_law_linear(ft, GF)
%FH_LAW_LINEAR Linear softening sigma = ft (1 - w/wc), wc = 2 GF/ft (Hillerborg 1976;
%   law used by Carpinteri, Cadamuro & Ventura 2015).
law.w = [0, 2*GF/ft]; law.s = [ft, 0]; law.tail = 'zero'; law.name = 'linear';
end
