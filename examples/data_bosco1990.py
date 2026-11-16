"""
Lightly reinforced concrete beams of Bosco, Carpinteri and Debernardi (1990), 2nd series.

Source: C. Bosco, A. Carpinteri, P.G. Debernardi, "Use of the brittleness number as a rational
approach to minimum reinforcement design", in: L. Elfgren, S.P. Shah (eds.), Analysis of Concrete
Structures by Fracture Mechanics, RILEM Proc. 6, Chapman and Hall, 1990, pp. 121-138.
Table 2 (geometry, span, bars, As/A, N_P, cracking/yielding/ultimate loads) and Table 3
(fcm, Ec, GF, KIC per concrete grade), transcribed without modification. Three-point bending.

Not reported in Table 2 and therefore DERIVED/ASSUMED here (declared in the documents):
  * fy: back-calculated from the reported N_P with Eq. (15): fy = N_P KIC / ((As/A) h^0.5);
  * concrete cover to the bar centroid: c = 0.1 h (value stated for the 1st series);
  * tensile strength: MC2010 fctm(fck = fcm - 8);
  * self-weight: 24 kN/m3 (included in the bending moment).
"""

GRADES = {  # Table 3
    "I": {"fcm": 16.63, "Ec": 21000.0, "GF": 0.143, "KIC": 54.80},
    "II": {"fcm": 29.38, "Ec": 23150.0, "GF": 0.134, "KIC": 55.70},
    "III": {"fcm": 31.96, "Ec": 23500.0, "GF": 0.137, "KIC": 56.74},
}

# id, grade, b, h, span, n_bars, phi, As/A [%], NP, P_crack, P_yield, P_ult [kN]
BEAMS = [
    ("A1", "I", 150, 100, 600, 1, 5, 0.131, 0.151, 5.95, 8.37, 8.95),
    ("A2", "I", 150, 100, 600, 2, 5, 0.261, 0.301, 5.87, 12.76, 16.29),
    ("A3", "I", 150, 100, 600, 3, 5, 0.392, 0.452, 5.76, 21.43, 24.04),
    ("AA1", "II", 150, 100, 600, 1, 5, 0.131, 0.148, 11.26, 11.49, 11.49),
    ("AA2", "II", 150, 100, 600, 2, 5, 0.261, 0.297, 11.10, 19.80, 21.21),
    ("AA3", "II", 150, 100, 600, 3, 5, 0.392, 0.445, 11.25, 27.40, 29.23),
    ("AAA1", "III", 200, 100, 600, 1, 5, 0.098, 0.110, 10.62, 8.75, 7.50),
    ("AAA2", "III", 200, 100, 600, 2, 5, 0.197, 0.220, 9.75, 17.50, 17.50),
    ("B1", "I", 150, 200, 1200, 1, 6, 0.094, 0.120, 9.67, 9.73, 10.99),
    ("B2", "I", 150, 200, 1200, 2, 6, 0.188, 0.237, 9.80, 15.69, 20.33),
    ("B3", "I", 150, 200, 1200, 3, 6, 0.283, 0.356, 9.30, 24.44, 28.89),
    ("BB1", "II", 150, 200, 1200, 1, 6, 0.094, 0.117, 11.34, 10.12, 12.26),
    ("BB2", "II", 150, 200, 1200, 2, 6, 0.188, 0.234, 14.05, 17.69, 23.35),
    ("BB3", "II", 150, 200, 1200, 3, 6, 0.283, 0.351, 12.11, 24.91, 33.00),
    ("BBB1", "III", 200, 200, 1200, 1, 6, 0.070, 0.086, 16.75, 9.75, 11.25),
    ("BBB2", "III", 200, 200, 1200, 2, 6, 0.141, 0.172, 18.95, 20.12, 24.70),
    ("C1", "I", 150, 400, 2400, 1, 8, 0.083, 0.145, 15.88, 17.69, 20.23),
    ("C2", "I", 150, 400, 2400, 2, 8, 0.167, 0.291, 18.46, 32.21, 37.20),
    ("C3", "I", 150, 400, 2400, 3, 8, 0.250, 0.436, 19.98, 47.59, 52.69),
    ("CC1", "II", 150, 400, 2400, 1, 8, 0.083, 0.144, 19.94, 20.06, 23.65),
    ("CC2", "II", 150, 400, 2400, 2, 8, 0.167, 0.287, 21.73, 31.31, 37.99),
    ("CC3", "II", 150, 400, 2400, 3, 8, 0.250, 0.431, 18.75, 47.50, 55.94),
    ("CCC1", "III", 150, 400, 2400, 1, 8, 0.083, 0.141, 18.25, 19.00, 23.10),
    ("CCC2", "III", 150, 400, 2400, 2, 8, 0.167, 0.283, 19.62, 28.75, 35.00),
    ("D2", "I", 200, 800, 4800, 2, 10, 0.098, 0.231, 30.05, 44.05, 51.48),
    ("D3", "I", 200, 800, 4800, 3, 10, 0.147, 0.346, 44.27, 63.92, 75.45),
    ("DD1", "II", 200, 800, 4800, 1, 10, 0.049, 0.113, 35.24, 27.78, 29.10),
    ("DD2", "II", 200, 800, 4800, 2, 10, 0.098, 0.227, 43.37, 45.33, 44.40),
    ("DD3", "II", 200, 800, 4800, 3, 10, 0.147, 0.340, 45.60, 60.60, 97.65),
]
