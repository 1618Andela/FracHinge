% EXAMPLE_CARPINTERI  R/SFRC beams of Ventura et al. reported by Carpinteri, Cadamuro & Ventura
% (2015, Mater. Struct. 48:235-247): b = 100, h = 200, L = 1200 mm, bars at c = 12 mm.
% Inputs B: ft and GF of Table 3, Ec = 22199 MPa, fy = fu (Table 8). Experimental peak loads
% digitised from Figs. 7-10 (mean of the curves). Compare with examples/generate_all.py.
b = 100; h = 200; L = 1200; c = 12; Ec = 22199; Es = 206000;
vf  = [0 10 20 40];
fc  = [35.88 50.95 56.85 45.90];
ft  = [2.8 3.7 4.0 3.4];
GF  = [0.113 0.441 0.492 1.772];
nb  = [0 1 2 2 2]; phi = [0 8 8 12 20]; fu = [0 630 630 646 625];
Pexp = [15.3 17.5 35.5 80.0 160.0; 16.0 21.5 40.0 85.0 150.0;
        17.0 25.0 44.0 88.0 170.0; 16.0 27.5 42.0 89.0 155.0];
Msw = 24e-6*b*h*L^2/8;
ratio = zeros(4, 5);
for i = 1:4
  conc = fh_concrete('ec2', fc(i), Ec);
  fctm = 0.3*(fc(i) - 8)^(2/3);
  for j = 1:5
    if nb(j) > 0
      rb = fh_rebar(c, nb(j)*pi*phi(j)^2/4, phi(j), 0.999*fu(j), fu(j), Es, 0.10, 2*fctm, fctm);
    else
      rb = [];
    end
    H = fh_hinge(b, h, fh_law_linear(ft(i), GF(i)), conc, Ec, h/2, rb, 300);
    R = fh_response(H, 0.06, 150);
    ok = logical(cumprod(-R.eps_top <= conc.eps_cu & R.y0 >= 0 & R.y0 <= h));
    P = 4*(max(R.M(ok)) - Msw)/L/1000;
    ratio(i, j) = P/Pexp(i, j);
    fprintf('Vf = %2d kg/m3, bars %d: P = %6.1f kN (exp %6.1f) ratio %.2f\n', vf(i), j, P, Pexp(i, j), ratio(i, j));
  end
end
fprintf('mean ratio %.2f, CoV %.1f %%\n', mean(ratio(:)), 100*std(ratio(:))/mean(ratio(:)));
