"""
Safety formats for non-linear (fracture-mechanics based) analyses — fib Model Code 2010, 7.11.3.

  * Partial factor method (7.11.3.4): design input values f_d -> R_d = r(f_d).
  * Global resistance factor (7.11.3.3.1): R_d = r(f_m*) / (gamma_R gamma_Rd), gamma_R = 1.2,
    gamma_Rd = 1.06, with reduced mean values (steel f_ym, concrete 0.85 f_cm).
  * ECOV (7.11.3.3.2): V_R = ln(R_m/R_k)/1.65, gamma_R = exp(alpha_R beta V_R),
    R_d = R_m / (gamma_R gamma_Rd), with alpha_R = 0.8 and beta = 3.8.
"""
from __future__ import annotations

import numpy as np


def ecov(R_m: float, R_k: float, alpha_R: float = 0.8, beta: float = 3.8, gamma_Rd: float = 1.06):
    VR = np.log(R_m / R_k) / 1.65
    gR = np.exp(alpha_R * beta * VR)
    return {"V_R": VR, "gamma_R": gR, "R_d": R_m / (gR * gamma_Rd), "gamma_Rd": gamma_Rd}


def global_resistance(R_mstar: float, gamma_R: float = 1.2, gamma_Rd: float = 1.06):
    return R_mstar / (gamma_R * gamma_Rd)
