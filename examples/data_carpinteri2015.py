"""R/SFRC beams tested by Ventura et al. (2003, Politecnico di Torino) and reported by
Carpinteri, Cadamuro & Ventura (2015), Materials and Structures 48:235-247.

b = 100 mm, h = 200 mm, length 1400 mm, three-point bending, bar layer at c = 12 mm
from the soffit (Sect. 3 of the paper). Span L = 1200 mm (Cortes Buitrago, 2023, Table 2.7).
Material data: Tables 3, 7 and 8 of the paper. Experimental peak loads were DIGITISED
from Figs. 7-10 of the paper (mean of the experimental curves; min-max range given);
uncertainty about +/- 1 kN for P < 50 kN and +/- 5 kN for P > 50 kN.
"""
B, H, L, C = 100.0, 200.0, 1200.0, 12.0
EC = 22199.0          # Table 4 / Table 8
ES = 206000.0
# fibre content: (fc [Table 3], ft Table 8 [notched flexural], GF Table 8, ft Table 3 [MC90], GF,tot Table 3)
MAT = {0: dict(fc=35.88, ft8=5.20, GF8=0.113, ft3=2.8, GF3=0.113),
       10: dict(fc=50.95, ft8=5.52, GF8=0.441, ft3=3.7, GF3=0.441),
       20: dict(fc=56.85, ft8=6.96, GF8=0.492, ft3=4.0, GF3=0.492),
       40: dict(fc=45.90, ft8=8.06, GF8=1.772, ft3=3.4, GF3=1.772)}
# bars: label -> (n, phi, fy, fu) [Table 7]
BARS = {"0Ø0": (0, 0.0, 0.0, 0.0), "1Ø8": (1, 8.0, 484.0, 630.0), "2Ø8": (2, 8.0, 484.0, 630.0),
        "2Ø12": (2, 12.0, 572.0, 646.0), "2Ø20": (2, 20.0, 503.0, 625.0)}
# digitised experimental peak loads [kN]: (mean, min, max)
PEXP = {
    0: {"0Ø0": (15.3, 13.0, 17.5), "1Ø8": (17.5, 16.5, 18.0), "2Ø8": (35.5, 34.0, 36.5),
        "2Ø12": (80.0, 76.0, 84.0), "2Ø20": (160.0, 150.0, 170.0)},
    10: {"0Ø0": (16.0, 12.0, 18.5), "1Ø8": (21.5, 20.0, 22.5), "2Ø8": (40.0, 38.0, 42.0),
         "2Ø12": (85.0, 80.0, 88.0), "2Ø20": (150.0, 140.0, 160.0)},
    20: {"0Ø0": (17.0, 15.0, 18.5), "1Ø8": (25.0, 23.0, 26.0), "2Ø8": (44.0, 42.0, 46.0),
         "2Ø12": (88.0, 80.0, 92.0), "2Ø20": (170.0, 160.0, 175.0)},
    40: {"0Ø0": (16.0, 14.0, 17.5), "1Ø8": (27.5, 26.0, 30.0), "2Ø8": (42.0, 40.0, 44.0),
         "2Ø12": (89.0, 85.0, 93.0), "2Ø20": (155.0, 130.0, 165.0)},
}
