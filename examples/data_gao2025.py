"""
Experimental data used for validation (transcribed from the published tables; no values invented).

Gao, D., et al. (2025). Flexural tensile property and damage evolution of hooked-end steel fiber
reinforced concrete with low fiber dosage and concrete strength. Journal of Building Engineering 99,
111680.  EN 14651 notched prisms 150x150x550 mm, notch 25 mm, span 500 mm; mean of 6 specimens.

Columns: id, cube compressive strength fcu [MPa] (Table 3), fR1..fR4 [MPa] (Table 4),
Fmax-1 [kN] (first peak, Table 5), E [GPa] (Table 5).
"""

GAO2025 = [
    # id,            fcu,   fR1,   fR2,   fR3,   fR4,   Fmax1,  E
    ("3D65/60-10", 42.57, 1.445, 1.128, 1.159, 1.123, 11.943, 26.973),
    ("3D65/60-20", 43.26, 2.155, 2.056, 2.083, 1.951, 12.032, 24.680),
    ("3D80/60-20", 42.50, 2.483, 2.746, 2.844, 2.772, 13.886, 22.168),
    ("3D100/60-20", 41.89, 2.726, 3.115, 3.482, 3.660, 14.826, 20.642),
    ("3D65/60-30", 45.50, 2.908, 3.089, 3.065, 2.930, 14.102, 24.247),
    ("3D65/60-40", 48.24, 3.651, 4.587, 4.817, 4.655, 14.610, 24.664),
    ("4D65/60-10", 43.37, 1.234, 1.222, 1.333, 1.268, 12.396, 21.028),
    ("4D65/60-20", 43.74, 1.967, 2.825, 3.187, 3.491, 13.560, 19.975),
    ("4D65/60-30", 46.37, 2.741, 3.675, 4.231, 4.520, 14.181, 19.715),
    ("4D65/60-40", 47.73, 4.697, 5.610, 6.094, 6.254, 14.658, 19.046),
    ("5D65/60-10", 42.96, 1.743, 1.973, 2.390, 2.560, 12.153, 22.894),
    ("5D65/60-20", 44.42, 2.162, 2.962, 3.500, 3.909, 15.030, 19.082),
    ("5D65/60-30", 46.08, 3.284, 4.651, 4.971, 5.004, 15.129, 22.034),
    ("5D65/60-40", 48.68, 4.949, 6.589, 6.857, 6.798, 16.345, 22.363),
]
# Plain matrix C30 (Table 3 and Table 5): fcu = 38.42 MPa, Fmax-1 = 11.827 kN, E = 23.895 GPa
GAO2025_PLAIN = ("C30", 38.42, 11.827, 23.895)
