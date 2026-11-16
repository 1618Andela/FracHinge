"""Slender RC and SFRC columns under constant axial load and cyclic lateral load (central stub,
semi-column length Ls = 1500 mm, b = 150 mm, 6 longitudinal bars, mechanical cover 0.15 h).
Sources:
  Caballero-Morrison, K.E., Bonet, J.L., Navarro-Gregori, J., Marti-Vargas, J.R. (2012). Behaviour of
    steel-fibre-reinforced normal-strength concrete slender columns under cyclic loading.
    Engineering Structures 39, 162-175 (Tables 1, 3, 5).
  Caballero-Morrison, K.E., Bonet, J.L., Navarro-Gregori, J., Serna-Ros, P. (2013). An experimental
    study of steel fiber-reinforced high-strength concrete slender columns under cyclic loading.
    Engineering Structures 57, 565-577 (Tables 1, 3, 4). Specimens N12-N15 of the 2013 paper repeat
    N1, N7, N9 and N11 of 2012 and are not duplicated here.
Columns: (source, id, h [mm], fc [MPa], N [kN], bar diameter [mm], stirrup spacing [mm],
          fibres [kg/m3], fR1 [MPa], fR3 [MPa], M_max [kN m] (critical section, second order and
          self weight included by the authors)).
Steel (both papers): phi12 fy = 548.27, fu = 656.98 MPa, eps_u = 0.184; phi8 fy = 545.52,
fu = 648.21 MPa, eps_u = 0.172; Es = 200 GPa.
"""
DB = [
    ("2012", "NF00L05V2S100", 260, 33.57, 491.7, 12, 100, 0, None, None, 71.68),
    ("2012", "NF00L05V2S50", 260, 29.68, 413.9, 12, 50, 0, None, None, 67.38),
    ("2012", "NF00L05V1S50", 260, 29.18, 138.2, 12, 50, 0, None, None, 47.14),
    ("2012", "NF00L10V2S70", 140, 30.63, 211.2, 8, 70, 0, None, None, 16.29),
    ("2012", "NF00L05V3S50", 260, 30.40, 637.9, 12, 50, 0, None, None, 69.55),
    ("2012", "NF00L10V3S70", 140, 32.68, 347.0, 8, 70, 0, None, None, 16.56),
    ("2012", "NF30L05V2S100", 260, 33.37, 473.8, 12, 100, 30, 2.89, 2.91, 71.25),
    ("2012", "NF30L05V2S50", 260, 38.00, 523.3, 12, 50, 30, 2.88, 3.48, 75.60),
    ("2012", "NF60L05V2S100", 260, 33.65, 412.7, 12, 100, 60, 6.51, 7.18, 70.83),
    ("2012", "NF60L05V2S50", 260, 29.52, 448.6, 12, 50, 60, 6.27, 6.77, 69.36),
    ("2012", "NF60L05V2S600", 260, 32.12, 420.6, 12, 600, 60, 6.59, 7.46, 71.24),
    ("2012", "NF60L10V2S70", 140, 30.35, 241.4, 8, 70, 60, 6.74, 7.52, 17.35),
    ("2012", "NF60L10V1S70", 140, 30.78, 107.9, 8, 70, 60, 7.79, 7.97, 14.83),
    ("2012", "NF60L05V1S50", 260, 32.43, 153.4, 12, 50, 60, 7.30, 8.15, 53.16),
    ("2013", "HF00L05V2S100", 260, 64.10, 895.1, 12, 100, 0, None, None, 114.63),
    ("2013", "HF00L05V2S50", 260, 71.70, 959.0, 12, 50, 0, None, None, 90.17),
    ("2013", "HF00L05V1S50", 260, 75.40, 294.7, 12, 50, 0, None, None, 65.59),
    ("2013", "HF00L10V1S70", 140, 70.50, 175.2, 8, 70, 0, None, None, 18.36),
    ("2013", "HF30L05V2S600", 260, 75.60, 967.7, 12, 600, 30, 5.20, 7.15, 118.27),
    ("2013", "HF60L05V2S100", 260, 75.50, 1013.3, 12, 100, 60, 9.59, 10.01, 123.39),
    ("2013", "HF60L05V2S50", 260, 79.20, 1056.7, 12, 50, 60, 11.29, 11.42, 111.47),
    ("2013", "HF60L05V2S600", 260, 75.10, 993.7, 12, 600, 60, 12.53, 12.52, 118.66),
    ("2013", "HF60L05V1S50", 260, 81.10, 333.9, 12, 50, 60, 11.05, 14.46, 78.37),
    ("2013", "HF60L10V1S70", 140, 79.30, 190.2, 8, 70, 60, 12.21, 13.61, 19.73),
    ("2013", "HF60L10V2S70", 140, 79.10, 598.8, 8, 70, 60, 9.41, 9.25, 32.41),
]
