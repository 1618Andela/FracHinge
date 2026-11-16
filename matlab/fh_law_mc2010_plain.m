function law = fh_law_mc2010_plain(fctm, GF)
%FH_LAW_MC2010_PLAIN Bilinear softening of plain concrete, fib MC2010 Eqs. (5.1-31)-(5.1-32).
%   law = fh_law_mc2010_plain(fctm [MPa], GF [N/mm]) -> struct with nodes w [mm], s [MPa].
w1 = GF/fctm; wc = 5*GF/fctm;
law.w = [0, w1, wc]; law.s = [fctm, 0.2*fctm, 0]; law.tail = 'zero';
law.name = 'MC2010 plain';
end
