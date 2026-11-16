function rb = fh_rebar(y, As, phi, fy, fu, Es, eps_u, tau0, tau1)
%FH_REBAR Bonded bar layer crossing the crack (Tension Chord Model bond: tau0 = 2 fctm,
%   tau1 = fctm; Marti et al. 1998).
rb.y = y; rb.As = As; rb.phi = phi; rb.fy = fy; rb.fu = fu; rb.Es = Es; rb.eps_u = eps_u;
rb.Esh = (fu - fy)/max(eps_u - fy/Es, 1e-9); rb.tau0 = tau0; rb.tau1 = tau1;
end
