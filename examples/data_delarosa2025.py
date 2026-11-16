"""De La Rosa, Ruiz & Carmona (2025), Engineering Fracture Mechanics 318:110918, analytical
planar-crack model for R/SFRC sections, validated on the beams of Carpinteri et al. (2015).
Table 1: material parameters used by the authors; Table 2: model moments (bond parameter eta = 4,
failure by w = wu) and experimental peak loads read by the authors (P_m,exp) with the experimental
moment M_m,exp = P_m,exp * 1.4 / 4 (i.e. computed with l = 1400 mm, the beam length)."""
# fibre content [kg/m3] -> (phi_f, fc, E [MPa], fR1, fR3, fFts, fFtu)
MAT = {10: (0.001, 51.0, 22200.0, 2.7, 3.2, 1.2, 1.1),
       20: (0.003, 56.9, 22200.0, 4.1, 4.6, 1.8, 1.5),
       40: (0.005, 45.9, 22200.0, 5.4, 5.6, 2.4, 1.7)}
# (fibre content, bars) -> (M_m model eta=4 [kNm], P_m,exp [kN], M_m,exp [kNm])
TAB2 = {(10, "1Ø8"): (6.2, 22.0, 7.7), (10, "2Ø8"): (10.2, 39.0, 13.7), (10, "2Ø12"): (20.0, 80.0, 28.0),
        (20, "1Ø8"): (7.0, 24.0, 8.4), (20, "2Ø8"): (11.0, 39.0, 13.7), (20, "2Ø12"): (20.6, 80.0, 28.0),
        (40, "1Ø8"): (7.6, 32.0, 11.2), (40, "2Ø8"): (11.5, 44.0, 15.4), (40, "2Ø12"): (21.2, 82.0, 28.7)}
