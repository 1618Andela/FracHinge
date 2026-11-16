"""Full-scale RC/SFRC beams of Minelli, Conforti, Cuenca & Plizzari (2014), Materials and
Structures 47:459-473, Tables 1, 4 and 5 (three-point loading, a/d = 3, b = 250 mm, d' = 60 mm,
8 bars in two layers of 4; fibres 50/0.8 hooked, 0 / 50 / 75 kg/m3).
"""
MAT = {"PC": dict(fcm=38.7, fctm=3.0, Ecm=33500.0, fR=None),
       "FRC50": dict(fcm=32.1, fctm=2.4, Ecm=30800.0, fR=(5.4, 5.6, 5.0, 4.5)),
       "FRC75": dict(fcm=33.1, fctm=2.5, Ecm=32100.0, fR=(6.0, 6.1, 6.0, 5.5))}
# H: (d, span, a, phi, As, fy, fu)
GEOM = {500: (440.0, 2640.0, 1320.0, 14.0, 1232.0, 506.0, 599.0),
        1000: (940.0, 5640.0, 2820.0, 20.0, 2513.0, 555.0, 651.0),
        1500: (1440.0, 8640.0, 4320.0, 24.0, 3619.0, 518.0, 612.0)}
# (H, mix): (failure, Vu kN, Mu kNm exp, Mu,fl kNm as computed by the authors)
TESTS = {(500, "PC"): ("S", 116, 153, 254), (500, "FRC50"): ("SY", 240, 316, 285),
         (500, "FRC75"): ("SY", 235, 310, 293), (1000, "PC"): ("S", 188, 529, 1210),
         (1000, "FRC50"): ("S", 272, 767, 1325), (1000, "FRC75"): ("S", 351, 989, 1356),
         (1500, "PC"): ("S", 211, 911, 2511), (1500, "FRC50"): ("S", 484, 2089, 2791),
         (1500, "FRC75"): ("S", 554, 2394, 2864)}
