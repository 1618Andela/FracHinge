"""
Fracture-mechanics quantities and size-effect tools.

  * Hillerborg et al. (1976): fictitious crack model; characteristic length lch = E GF / ft^2.
  * Carpinteri (1984), Bosco & Carpinteri (1992): brittleness numbers s_E = KIC/(ft h^0.5) and
    N_P = fy rho h^0.5 / KIC for reinforced members.
  * Bazant (1984): size effect law sigma_N = B ft / sqrt(1 + D/D0).
  * Tada, Paris & Irwin (2000) / Srawley (1976): K_I of a single-edge notched beam in 3PB, S/W = 4.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

from .materials import characteristic_length


def kic_from_gf(GF: float, E: float) -> float:
    """Irwin: K_IC = sqrt(E G_F) (plane stress) [MPa sqrt(mm)]."""
    return np.sqrt(E * GF)


def brittleness_hillerborg(h: float, E: float, GF: float, ft: float) -> float:
    """Hillerborg's brittleness number h / lch."""
    return h / characteristic_length(E, GF, ft)


def brittleness_carpinteri_sE(KIC: float, ft: float, h: float) -> float:
    """Carpinteri's energy brittleness number s_E = KIC / (ft h^0.5)."""
    return KIC / (ft * np.sqrt(h))


def brittleness_NP(fy: float, rho: float, h: float, KIC: float) -> float:
    """Bosco-Carpinteri reinforcement brittleness number N_P = fy rho h^0.5 / KIC."""
    return fy * rho * np.sqrt(h) / KIC


def sel_bazant(D, B_ft, D0):
    """Bazant's size-effect law sigma_N = B ft / sqrt(1 + D/D0)."""
    return B_ft / np.sqrt(1.0 + np.asarray(D) / D0)


def fit_sel(D, sigmaN):
    p, _ = curve_fit(sel_bazant, np.asarray(D, float), np.asarray(sigmaN, float),
                     p0=[max(sigmaN) * 1.2, np.median(D)], maxfev=20000)
    return {"B_ft": p[0], "D0": p[1]}


def senb_geometry_factor(alpha: float) -> float:
    """F(alpha) for SENB, S/W = 4: K_I = sigma_N sqrt(a) F(alpha), sigma_N = 6M/(B W^2)."""
    a = alpha
    return (1.99 - a * (1 - a) * (2.15 - 3.93 * a + 2.7 * a ** 2)) / ((1 + 2 * a) * (1 - a) ** 1.5)


def k1_senb(M: float, B: float, W: float, a: float) -> float:
    sigma = 6.0 * M / (B * W ** 2)
    return sigma * np.sqrt(a) * senb_geometry_factor(a / W)
